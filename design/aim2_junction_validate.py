#!/usr/bin/env python3
"""
Aim 2 – ProteinMPNN junction optimization validation.

Pipeline:
  1. Parse 20 MPNN junction sequences from mpnn_junctions/mpnn_output/seqs/complex.fa
  2. Run Boltz-2 FULL (5-chain) for any not already done
  3. Score BRD4 contact fractions (frag residues only) from FULL PDBs
  4. Pick top 5 by FULL BRD4 CF
  5. Run Boltz-2 NO_A + NO_B for those top 5
  6. Report AND-gate scores; save TSV

AND-gate criterion: FULL CF > 0.10  and  NO_A CF < 0.05  and  NO_B CF < 0.05

Usage:
  python design/aim2_junction_validate.py [--top-n 5] [--status]
"""
import argparse, json, subprocess, csv
from pathlib import Path
import numpy as np

PROJECT   = Path("/mnt/d/gbm-ANDGATE")
BOLTZ     = Path("/home/ayushregmibagale/miniconda3/envs/boltz2/bin/boltz")
MPNN_FA   = PROJECT / "design" / "mpnn_junctions" / "mpnn_output" / "seqs" / "complex.fa"
OUT_DIR   = PROJECT / "design" / "aim2_boltz"
RESULTS   = PROJECT / "design" / "aim2_boltz" / "aim2_results.tsv"

SOX2_HMG  = "DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLM"
SOX9_HMG  = "PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"
BRD4_BD1  = "NPPPPETSNPNKPKRQTNQLQYLLRVVLKTLWKHQFAWPFQQPVDAVKLNLPDYYKIIKTPMDMGTIKKRLENNYYWNAQECIQDFNTMFTNCYIYNKPGDDIVLMAEALEKLFLQKINELPTEE"
CONTACT_CUTOFF = 8.0

# frag_N = CompA residues 91-109 (19 aa BRD4 frag); frag_C = CompB residues 87-143 (57 aa)
FRAGNSEQ_START, FRAGNSEQ_END = 91, 109
FRAGCSEQ_START, FRAGCSEQ_END = 87, 143


def parse_mpnn_seqs(fa):
    """Return list of (id, comp_a, comp_b) from the MPNN FASTA (/ separates chains)."""
    seqs = []
    cur_id, cur_seq = None, None
    for line in open(fa):
        line = line.rstrip()
        if line.startswith(">"):
            if cur_id and cur_seq:
                seqs.append((cur_id, cur_seq))
            cur_id = line[1:].split(",")[0].strip()
            cur_seq = ""
        else:
            cur_seq += line
    if cur_id and cur_seq:
        seqs.append((cur_id, cur_seq))
    result = []
    for sid, seq in seqs:
        parts = seq.split("/")
        if len(parts) < 2:
            continue
        comp_a, comp_b = parts[0], parts[1]
        result.append((sid, comp_a, comp_b))
    # First entry is the native sequence (reference); skip it; samples start at index 1
    samples = result[1:]
    return samples


def write_yaml(path, chains):
    blocks = [f"  - protein:\n      id: {cid}\n      sequence: {seq}\n      msa: empty"
              for cid, seq in chains.items()]
    path.write_text("sequences:\n" + "\n".join(blocks) + "\n")


def _pred_dir(work_dir):
    """Boltz-2 writes to <out_dir>/boltz_results_input/predictions/input/."""
    nested = work_dir / "boltz_results_input" / "predictions" / "input"
    flat   = work_dir / "predictions" / "input"
    return nested if nested.exists() else flat


def run_boltz(work_dir, yaml_path):
    pred_dir = _pred_dir(work_dir)
    if pred_dir.exists() and any(pred_dir.glob("*.pdb")):
        return
    work_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [str(BOLTZ), "predict", str(yaml_path),
         "--out_dir", str(work_dir),
         "--recycling_steps", "3", "--sampling_steps", "200",
         "--output_format", "pdb", "--no_kernels"],
        capture_output=True, timeout=1200,
    )


def get_iptm(work_dir):
    cf = _pred_dir(work_dir) / "confidence_input_model_0.json"
    if not cf.exists():
        return 0.0
    d = json.load(open(cf))
    return float(d.get("iptm", 0.0))


def get_pdb(work_dir):
    pdbs = list(_pred_dir(work_dir).glob("input_model_0.pdb"))
    return str(pdbs[0]) if pdbs else ""


def parse_ca(pdb_path, chain, res_start=None, res_end=None):
    coords = {}
    try:
        with open(pdb_path) as f:
            for line in f:
                if not (line.startswith("ATOM") and line[12:16].strip() == "CA"):
                    continue
                if line[21] != chain:
                    continue
                resnum = int(line[22:26])
                if res_start is not None and resnum < res_start:
                    continue
                if res_end is not None and resnum > res_end:
                    continue
                coords[resnum] = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    except Exception:
        pass
    return coords


def cf_against(query_coords, target_coords, cutoff=CONTACT_CUTOFF):
    if not query_coords or not target_coords:
        return 0.0
    tgt = np.array(list(target_coords.values()))
    contacts = sum(
        np.linalg.norm(tgt - xyz, axis=1).min() < cutoff
        for xyz in query_coords.values()
    )
    return contacts / len(query_coords)


def compute_full_brd4_cf(pdb_path):
    """Combined frag_N+frag_C vs BRD4 CF in FULL 5-chain PDB."""
    # Boltz-2 FULL chain order: A=CompA, B=CompB, C=SOX2, D=SOX9, E=BRD4
    ca_fragn  = parse_ca(pdb_path, 'A', FRAGNSEQ_START, FRAGNSEQ_END)
    ca_fragc  = parse_ca(pdb_path, 'B', FRAGCSEQ_START, FRAGCSEQ_END)
    ca_brd4   = parse_ca(pdb_path, 'E')
    if not ca_brd4:
        return 0.0
    tgt = np.array(list(ca_brd4.values()))
    all_frags = list(ca_fragn.values()) + list(ca_fragc.values())
    if not all_frags:
        return 0.0
    contacts = sum(np.linalg.norm(tgt - xyz, axis=1).min() < CONTACT_CUTOFF for xyz in all_frags)
    return contacts / len(all_frags)


def compute_nob_brd4_cf(pdb_path):
    """frag_N vs BRD4 CF in NO_B 4-chain PDB (CompA only).
    Boltz-2 NO_B chain order: A=CompA, B=SOX2, C=SOX9, D=BRD4"""
    ca_fragn = parse_ca(pdb_path, 'A', FRAGNSEQ_START, FRAGNSEQ_END)
    ca_brd4  = parse_ca(pdb_path, 'D')
    return cf_against(ca_fragn, ca_brd4)


def compute_noa_brd4_cf(pdb_path):
    """frag_C vs BRD4 CF in NO_A 4-chain PDB (CompB only).
    Boltz-2 NO_A chain order: A=CompB, B=SOX2, C=SOX9, D=BRD4"""
    ca_fragc = parse_ca(pdb_path, 'A', FRAGCSEQ_START, FRAGCSEQ_END)
    ca_brd4  = parse_ca(pdb_path, 'D')
    return cf_against(ca_fragc, ca_brd4)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-n", type=int, default=5, help="Number of top seqs for NO_A/NO_B validation")
    parser.add_argument("--status", action="store_true", help="Show progress without running")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    samples = parse_mpnn_seqs(MPNN_FA)
    print(f"Parsed {len(samples)} MPNN junction sequences")

    # --- Stage 1: FULL Boltz-2 for all samples ---
    full_results = []
    for idx, (sid, comp_a, comp_b) in enumerate(samples, 1):
        tag = f"seq{idx:02d}"
        wd  = OUT_DIR / tag / "full"
        wd.mkdir(parents=True, exist_ok=True)

        # Also check old boltz_seq dirs (seq01-10 already run there)
        old_wd = PROJECT / "design" / "mpnn_junctions" / f"boltz_seq{idx:02d}" / "boltz_results_input"
        if old_wd.exists() and any(old_wd.glob("**/*.pdb")):
            # Use old results; copy confidence json location
            old_pred = old_wd / "predictions" / "input"
            iptm = get_iptm(old_wd)
            pdbs = list(old_pred.glob("*.pdb"))
            pdb  = str(pdbs[0]) if pdbs else ""
        else:
            if not args.status:
                yaml = wd / "input.yaml"
                write_yaml(yaml, {"A": comp_a, "B": comp_b,
                                  "C": SOX2_HMG, "D": SOX9_HMG, "E": BRD4_BD1})
                print(f"  FULL {tag} ...", flush=True)
                run_boltz(wd, yaml)
            iptm = get_iptm(wd)
            pdb  = get_pdb(wd)

        cf_full = compute_full_brd4_cf(pdb) if pdb else 0.0
        full_results.append({
            "id": tag, "comp_a": comp_a, "comp_b": comp_b,
            "iptm_full": iptm, "cf_brd4_full": cf_full, "pdb_full": pdb,
        })
        print(f"  {tag}: ipTM={iptm:.3f}  BRD4-CF(full)={cf_full:.3f}")

    # --- Stage 2: Top-N by FULL BRD4 CF → NO_A + NO_B ---
    full_results.sort(key=lambda r: r["cf_brd4_full"], reverse=True)
    top = full_results[:args.top_n]
    print(f"\nTop {args.top_n} by FULL BRD4 CF:")
    for r in top:
        print(f"  {r['id']}: BRD4-CF={r['cf_brd4_full']:.3f}  ipTM={r['iptm_full']:.3f}")

    all_rows = []
    for r in full_results:
        row = {**r, "iptm_nob": 0.0, "cf_brd4_nob": 0.0,
               "iptm_noa": 0.0, "cf_brd4_noa": 0.0, "pass": False}
        all_rows.append(row)

    for r in top:
        tag = r["id"]
        comp_a, comp_b = r["comp_a"], r["comp_b"]

        # NO_B (CompA only)
        wd_nob = OUT_DIR / tag / "no_b"
        wd_nob.mkdir(parents=True, exist_ok=True)
        if not args.status:
            yaml = wd_nob / "input.yaml"
            write_yaml(yaml, {"A": comp_a, "B": SOX2_HMG, "C": SOX9_HMG, "D": BRD4_BD1})
            print(f"  NO_B {tag} ...", flush=True)
            run_boltz(wd_nob, yaml)
        iptm_nob = get_iptm(wd_nob)
        pdb_nob  = get_pdb(wd_nob)
        cf_nob   = compute_nob_brd4_cf(pdb_nob) if pdb_nob else 0.0

        # NO_A (CompB only)
        wd_noa = OUT_DIR / tag / "no_a"
        wd_noa.mkdir(parents=True, exist_ok=True)
        if not args.status:
            yaml = wd_noa / "input.yaml"
            write_yaml(yaml, {"A": comp_b, "B": SOX2_HMG, "C": SOX9_HMG, "D": BRD4_BD1})
            print(f"  NO_A {tag} ...", flush=True)
            run_boltz(wd_noa, yaml)
        iptm_noa = get_iptm(wd_noa)
        pdb_noa  = get_pdb(wd_noa)
        cf_noa   = compute_noa_brd4_cf(pdb_noa) if pdb_noa else 0.0

        passed = (r["cf_brd4_full"] > 0.10 and cf_nob < 0.05 and cf_noa < 0.05)
        print(f"  {tag}: FULL={r['cf_brd4_full']:.3f}  NO_B={cf_nob:.3f}  NO_A={cf_noa:.3f}  "
              f"{'PASS' if passed else 'fail'}")

        for row in all_rows:
            if row["id"] == tag:
                row.update(iptm_nob=iptm_nob, cf_brd4_nob=cf_nob,
                           iptm_noa=iptm_noa, cf_brd4_noa=cf_noa, pass_=passed)
                break

    # --- Summary ---
    print("\n=== Aim 2 Junction Summary ===")
    print(f"{'ID':>6}  {'FULL-CF':>7}  {'NO_B-CF':>7}  {'NO_A-CF':>7}  {'FULL-ipTM':>9}  {'AND-gate':>8}")
    print("-" * 60)
    for r in full_results:
        nob = r.get("cf_brd4_nob", 0.0)
        noa = r.get("cf_brd4_noa", 0.0)
        passed = (r["cf_brd4_full"] > 0.10 and nob < 0.05 and noa < 0.05)
        flag = " ✓" if passed else ""
        print(f"{r['id']:>6}  {r['cf_brd4_full']:>7.3f}  {nob:>7.3f}  {noa:>7.3f}  "
              f"{r['iptm_full']:>9.3f}  {flag}")

    # Save TSV
    fields = ["id", "iptm_full", "cf_brd4_full", "iptm_nob", "cf_brd4_nob",
              "iptm_noa", "cf_brd4_noa", "comp_a", "comp_b"]
    with open(RESULTS, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rows)
    print(f"\nSaved → {RESULTS}")


if __name__ == "__main__":
    main()

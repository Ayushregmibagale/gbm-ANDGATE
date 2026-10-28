#!/usr/bin/env python3
"""
Aim 2b – Constrained ProteinMPNN junction optimization.

Fixes ALL BRD4-fragment residues; only redesigns GGS3 linker + binder C-terminus.

CompA chain layout (109 aa):
  NLS 1-10 | SOX2_binder 11-81 | GGS3 82-90 | frag_N 91-109
  Fixed:   1-77, 91-109   →   Redesigned: 78-90 (13 residues)

CompB chain layout (143 aa):
  NLS 1-10 | SOX9_binder 11-77 | GGS3 78-86 | frag_C 87-143
  Fixed:   1-74, 87-143   →   Redesigned: 75-86 (12 residues)

Chains C (SOX2), D (SOX9), E (BRD4): all fixed.

Usage:
  python design/aim2b_constrained_mpnn.py        # run MPNN + full validation
  python design/aim2b_constrained_mpnn.py --status
"""
import argparse, json, subprocess, csv, sys
from pathlib import Path
import numpy as np

PROJECT    = Path("/mnt/d/gbm-ANDGATE")
BOLTZ      = Path("/home/ayushregmibagale/miniconda3/envs/boltz2/bin/boltz")
COMPLEX_PDB = PROJECT / "design" / "mpnn_junctions" / "mpnn_output" / "complex.pdb"
OUT_DIR    = PROJECT / "design" / "aim2b"
MPNN_OUT   = OUT_DIR / "mpnn_output"
RESULTS    = OUT_DIR / "aim2b_results.tsv"

SOX2_HMG   = "DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLM"
SOX9_HMG   = "PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"
BRD4_BD1   = "NPPPPETSNPNKPKRQTNQLQYLLRVVLKTLWKHQFAWPFQQPVDAVKLNLPDYYKIIKTPMDMGTIKKRLENNYYWNAQECIQDFNTMFTNCYIYNKPGDDIVLMAEALEKLFLQKINELPTEE"

CONTACT_CUTOFF = 8.0
FRAGNSEQ_START, FRAGNSEQ_END   = 91, 109
FRAGCSEQ_START, FRAGCSEQ_END   = 87, 143


# ── ProteinMPNN ─────────────────────────────────────────────────────────────
def run_mpnn(n_seq=20, temperature=0.1, seed=42):
    MPNN_OUT.mkdir(parents=True, exist_ok=True)
    fa = MPNN_OUT / "seqs" / "complex.fa"
    if fa.exists() and fa.stat().st_size > 0:
        print(f"  MPNN already done ({fa})")
        return

    MPNN_PYTHON = Path("/home/ayushregmibagale/miniconda3/envs/proteinmpnn/bin/python")
    MPNN_PY_DIR = Path("/home/ayushregmibagale/projects/ProteinMPNN")

    parsed   = MPNN_OUT / "parsed_pdbs.jsonl"
    assigned = MPNN_OUT / "assigned_chains.jsonl"
    fixed    = MPNN_OUT / "fixed_positions.jsonl"

    # Parse PDB (helper needs a folder)
    pdb_dir = MPNN_OUT / "pdb_input"
    pdb_dir.mkdir(exist_ok=True)
    link = pdb_dir / "complex.pdb"
    if not link.exists():
        link.symlink_to(COMPLEX_PDB.resolve())

    def mpnn_run(script, extra_args):
        r = subprocess.run(
            [str(MPNN_PYTHON), str(MPNN_PY_DIR / "helper_scripts" / script)] + extra_args,
            capture_output=False)
        if r.returncode != 0:
            print(f"  {script} failed (rc={r.returncode})")
            sys.exit(1)

    print("  Parsing PDB ...", flush=True)
    mpnn_run("parse_multiple_chains.py",
             ["--input_path", str(pdb_dir), "--output_path", str(parsed)])

    print("  Assigning designable chains A B ...", flush=True)
    mpnn_run("assign_fixed_chains.py",
             ["--input_path", str(parsed), "--output_path", str(assigned),
              "--chain_list", "A B"])

    # Specify residues to REDESIGN (not fix) via --specify_non_fixed
    # A 78-90: GGS3 + last 4 of SOX2-binder (fixes all frag_N 91-109)
    # B 75-86: GGS3 + last 3 of SOX9-binder (fixes all frag_C 87-143)
    redesign_A = " ".join(str(i) for i in range(78, 91))
    redesign_B = " ".join(str(i) for i in range(75, 87))
    print("  Generating fixed positions (frag residues locked) ...", flush=True)
    mpnn_run("make_fixed_positions_dict.py",
             ["--input_path", str(parsed), "--output_path", str(fixed),
              "--chain_list", "A B",
              "--position_list", f"{redesign_A}, {redesign_B}",
              "--specify_non_fixed"])

    print("  Running ProteinMPNN ...", flush=True)
    r = subprocess.run(
        [str(MPNN_PYTHON), str(MPNN_PY_DIR / "protein_mpnn_run.py"),
         "--jsonl_path",            str(parsed),
         "--chain_id_jsonl",        str(assigned),
         "--fixed_positions_jsonl", str(fixed),
         "--out_folder",            str(MPNN_OUT),
         "--num_seq_per_target",    str(n_seq),
         "--sampling_temp",         str(temperature),
         "--seed",                  str(seed),
         "--batch_size",            "1"],
        capture_output=False)
    if r.returncode != 0:
        print(f"  MPNN failed (rc={r.returncode})")
        sys.exit(1)


def parse_mpnn_seqs(fa):
    seqs, cur_id, cur_seq = [], None, None
    for line in open(fa):
        line = line.rstrip()
        if line.startswith(">"):
            if cur_id and cur_seq:
                seqs.append((cur_id, cur_seq))
            cur_id, cur_seq = line[1:].split(",")[0].strip(), ""
        else:
            cur_seq += line
    if cur_id and cur_seq:
        seqs.append((cur_id, cur_seq))
    result = []
    for sid, seq in seqs:
        parts = seq.split("/")
        if len(parts) >= 2:
            result.append((sid, parts[0], parts[1]))
    return result[1:]   # skip native reference


# ── Boltz-2 helpers ─────────────────────────────────────────────────────────
def _pred_dir(work_dir):
    nested = work_dir / "boltz_results_input" / "predictions" / "input"
    flat   = work_dir / "predictions" / "input"
    return nested if nested.exists() else flat


def write_yaml(path, chains):
    blocks = [f"  - protein:\n      id: {cid}\n      sequence: {seq}\n      msa: empty"
              for cid, seq in chains.items()]
    path.write_text("sequences:\n" + "\n".join(blocks) + "\n")


def run_boltz(work_dir, yaml_path):
    pd = _pred_dir(work_dir)
    if pd.exists() and any(pd.glob("*.pdb")):
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
    return float(json.load(open(cf)).get("iptm", 0.0))


def get_pdb(work_dir):
    pdbs = list(_pred_dir(work_dir).glob("input_model_0.pdb"))
    return str(pdbs[0]) if pdbs else ""


def parse_ca(pdb_path, chain, r0=None, r1=None):
    coords = {}
    try:
        with open(pdb_path) as f:
            for line in f:
                if not (line.startswith("ATOM") and line[12:16].strip() == "CA"):
                    continue
                if line[21] != chain:
                    continue
                rn = int(line[22:26])
                if r0 is not None and rn < r0: continue
                if r1 is not None and rn > r1: continue
                coords[rn] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    except Exception:
        pass
    return coords


def cf_against(q_ca, t_ca, cutoff=CONTACT_CUTOFF):
    if not q_ca or not t_ca:
        return 0.0
    tgt = np.array(list(t_ca.values()))
    return sum(np.linalg.norm(tgt - xyz, axis=1).min() < cutoff for xyz in q_ca.values()) / len(q_ca)


def compute_full_brd4_cf(pdb):
    ca_fn = parse_ca(pdb, 'A', FRAGNSEQ_START, FRAGNSEQ_END)
    ca_fc = parse_ca(pdb, 'B', FRAGCSEQ_START, FRAGCSEQ_END)
    ca_b4 = parse_ca(pdb, 'E')
    if not ca_b4: return 0.0
    tgt = np.array(list(ca_b4.values()))
    frags = list(ca_fn.values()) + list(ca_fc.values())
    if not frags: return 0.0
    return sum(np.linalg.norm(tgt - xyz, axis=1).min() < CONTACT_CUTOFF for xyz in frags) / len(frags)


def compute_nob_brd4_cf(pdb):   # NO_B: A=CompA, D=BRD4
    return cf_against(parse_ca(pdb, 'A', FRAGNSEQ_START, FRAGNSEQ_END), parse_ca(pdb, 'D'))


def compute_noa_brd4_cf(pdb):   # NO_A: A=CompB, D=BRD4
    return cf_against(parse_ca(pdb, 'A', FRAGCSEQ_START, FRAGCSEQ_END), parse_ca(pdb, 'D'))


# ── main ────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--status", action="store_true")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Run MPNN
    if not args.status:
        run_mpnn(n_seq=20, temperature=0.1, seed=42)

    fa = MPNN_OUT / "seqs" / "complex.fa"
    if not fa.exists():
        print("MPNN output not found — run without --status first")
        return
    samples = parse_mpnn_seqs(fa)
    print(f"Parsed {len(samples)} constrained MPNN sequences")

    # 2. FULL Boltz-2 for all
    full_results = []
    for idx, (sid, comp_a, comp_b) in enumerate(samples, 1):
        tag = f"seq{idx:02d}"
        wd  = OUT_DIR / tag / "full"
        wd.mkdir(parents=True, exist_ok=True)
        if not args.status:
            yaml = wd / "input.yaml"
            write_yaml(yaml, {"A": comp_a, "B": comp_b,
                               "C": SOX2_HMG, "D": SOX9_HMG, "E": BRD4_BD1})
            print(f"  FULL {tag} ...", flush=True)
            run_boltz(wd, yaml)
        iptm = get_iptm(wd)
        pdb  = get_pdb(wd)
        cf   = compute_full_brd4_cf(pdb) if pdb else 0.0
        full_results.append({"id": tag, "comp_a": comp_a, "comp_b": comp_b,
                              "iptm_full": iptm, "cf_brd4_full": cf})
        print(f"  {tag}: ipTM={iptm:.3f}  BRD4-CF(full)={cf:.3f}")

    # 3. Top-N NO_A + NO_B
    full_results.sort(key=lambda r: r["cf_brd4_full"], reverse=True)
    top = full_results[:args.top_n]
    print(f"\nTop {args.top_n} by FULL BRD4 CF:")
    for r in top:
        print(f"  {r['id']}: BRD4-CF={r['cf_brd4_full']:.3f}  ipTM={r['iptm_full']:.3f}")

    all_rows = [{**r, "iptm_nob": 0.0, "cf_brd4_nob": 0.0,
                 "iptm_noa": 0.0, "cf_brd4_noa": 0.0, "and_gate": False}
                for r in full_results]

    for r in top:
        tag, comp_a, comp_b = r["id"], r["comp_a"], r["comp_b"]

        wd_nob = OUT_DIR / tag / "no_b"
        wd_nob.mkdir(parents=True, exist_ok=True)
        if not args.status:
            yaml = wd_nob / "input.yaml"
            write_yaml(yaml, {"A": comp_a, "B": SOX2_HMG, "C": SOX9_HMG, "D": BRD4_BD1})
            print(f"  NO_B {tag} ...", flush=True)
            run_boltz(wd_nob, yaml)
        cf_nob = compute_nob_brd4_cf(get_pdb(wd_nob)) if get_pdb(wd_nob) else 0.0

        wd_noa = OUT_DIR / tag / "no_a"
        wd_noa.mkdir(parents=True, exist_ok=True)
        if not args.status:
            yaml = wd_noa / "input.yaml"
            write_yaml(yaml, {"A": comp_b, "B": SOX2_HMG, "C": SOX9_HMG, "D": BRD4_BD1})
            print(f"  NO_A {tag} ...", flush=True)
            run_boltz(wd_noa, yaml)
        cf_noa = compute_noa_brd4_cf(get_pdb(wd_noa)) if get_pdb(wd_noa) else 0.0

        passed = r["cf_brd4_full"] > 0.10 and cf_nob < 0.05 and cf_noa < 0.05
        print(f"  {tag}: FULL={r['cf_brd4_full']:.3f}  NO_B={cf_nob:.3f}  NO_A={cf_noa:.3f}  "
              f"{'PASS ✓' if passed else 'fail'}")
        for row in all_rows:
            if row["id"] == tag:
                row.update(iptm_nob=get_iptm(wd_nob), cf_brd4_nob=cf_nob,
                           iptm_noa=get_iptm(wd_noa), cf_brd4_noa=cf_noa, and_gate=passed)
                break

    # 4. Summary
    print("\n=== Aim 2b Constrained MPNN Summary ===")
    print(f"{'ID':>6}  {'FULL-CF':>7}  {'NO_B-CF':>7}  {'NO_A-CF':>7}  {'ipTM':>6}  {'AND-gate':>8}")
    print("-" * 58)
    for r in all_rows:
        flag = " ✓" if r.get("and_gate") else ""
        print(f"{r['id']:>6}  {r['cf_brd4_full']:>7.3f}  {r.get('cf_brd4_nob',0):>7.3f}  "
              f"{r.get('cf_brd4_noa',0):>7.3f}  {r['iptm_full']:>6.3f}  {flag}")

    fields = ["id", "iptm_full", "cf_brd4_full", "iptm_nob", "cf_brd4_nob",
              "iptm_noa", "cf_brd4_noa", "and_gate", "comp_a", "comp_b"]
    with open(RESULTS, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rows)
    print(f"\nSaved → {RESULTS}")


if __name__ == "__main__":
    main()

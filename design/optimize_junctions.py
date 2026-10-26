#!/usr/bin/env python3
"""
ProteinMPNN junction optimization for AND-gate components.

Uses the Boltz-2 predicted 5-chain FULL complex structure to redesign
the GGS3 linker + 3 flanking residues at each SOX_binder–BRD4_frag junction.

Chain map in the Boltz-2 FULL PDB:
  A = Component A (109 aa): NLS(1-10) + SOX2binder(11-81) + GGS3(82-90) + BRD4fragN(91-109)
  B = Component B (143 aa): NLS(1-10) + SOX9binder(11-77) + GGS3(78-86) + BRD4fragC(87-143)
  C = SOX2_HMG  (fixed target)
  D = SOX9_HMG  (fixed target)
  E = BRD4_BD1  (fixed target)

Designable residues (everything else fixed):
  Chain A: 79-93  (GGS3 82-90 + 3 flanking each side)
  Chain B: 75-89  (GGS3 78-86 + 3 flanking each side)

Produces 20 redesigned sequence pairs, validated with Boltz-2 (FULL 5-chain).
"""

import csv
import json
import subprocess
import shutil
from pathlib import Path
import numpy as np

PROJECT  = Path("/mnt/d/gbm-ANDGATE")
RIPTAC   = Path("/mnt/d/gbm-RIPTAC")
MPNN_SCRIPT = Path("/home/ayushregmibagale/projects/ProteinMPNN/protein_mpnn_run.py")
MPNN_PYTHON = Path("/home/ayushregmibagale/miniconda3/envs/proteinmpnn/bin/python")
BOLTZ       = Path("/home/ayushregmibagale/miniconda3/envs/boltz2/bin/boltz")

FULL_PDB = RIPTAC / "and_gate/boltz_outputs/v1_FULL/boltz_results_input/predictions/input/input_model_0.pdb"
OUT_DIR  = PROJECT / "design" / "mpnn_junctions"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Original sequences
NLS         = "PKKKRKVGGS"
GGS3        = "GGSGGSGGS"
SOX2_BINDER = "PYEIPEEDKEWDEYSFAMVEVHERLMEAERQLRLLSPEEEEKKRELEEEIERARKELEEVHKRGEAIERKM"
SOX9_BINDER = "MEEVRELLERARELEEKVRELEDDPSKREERERLLEEAERLREEALDLLVEQIREPLLQWMREEGLL"
BRD4_FRAG_N = "MNKSIEELLEELRKKFKEN"
BRD4_FRAG_C = "TEENKEKMLKELEERSEQYKEAALKFGESEEKAEEVKEYSLKLNEQVKEGWKANNQK"

SOX2_HMG = "DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLM"
SOX9_HMG = "PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"
BRD4_BD1 = "NPPPPETSNPNKPKRQTNQLQYLLRVVLKTLWKHQFAWPFQQPVDAVKLNLPDYYKIIKTPMDMGTIKKRLENNYYWNAQECIQDFNTMFTNCYIYNKPGDDIVLMAEALEKLFLQKINELPTEE"

N_SEQS      = 20   # sequences per design run
BOLTZ_TOP_N = 10   # validate top N by MPNN score
IPTM_THRESH = 0.30
CF_THRESH   = 0.08

CONTACT_CUTOFF = 8.0


# ── Helpers ───────────────────────────────────────────────────────────────────

def parse_ca(pdb_path, chain):
    coords = {}
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip()=="CA" and line[21]==chain:
                    coords[int(line[22:26])] = np.array(
                        [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    except Exception:
        pass
    return coords


def brd4_cf(pdb_path, comp_chain, brd4_chain):
    ca_c = parse_ca(pdb_path, comp_chain)
    ca_b = parse_ca(pdb_path, brd4_chain)
    if not ca_c or not ca_b:
        return 0.0
    brd4 = np.array(list(ca_b.values()))
    return sum(np.linalg.norm(brd4 - xyz, axis=1).min() < CONTACT_CUTOFF
               for xyz in ca_c.values()) / len(ca_c)


def run_boltz_5chain(comp_a_seq, comp_b_seq, work_dir):
    work_dir.mkdir(parents=True, exist_ok=True)
    yaml = work_dir / "input.yaml"
    yaml.write_text(f"""sequences:
  - protein:
      id: A
      sequence: {comp_a_seq}
      msa: empty
  - protein:
      id: B
      sequence: {comp_b_seq}
      msa: empty
  - protein:
      id: C
      sequence: {SOX2_HMG}
      msa: empty
  - protein:
      id: D
      sequence: {SOX9_HMG}
      msa: empty
  - protein:
      id: E
      sequence: {BRD4_BD1}
      msa: empty
""")
    pred = work_dir / "boltz_results_input" / "predictions" / "input"
    if not pred.exists():
        subprocess.run([str(BOLTZ), "predict", str(yaml),
                        "--out_dir", str(work_dir),
                        "--recycling_steps", "3", "--sampling_steps", "200",
                        "--output_format", "pdb", "--no_kernels"],
                       capture_output=True, timeout=1200)
    sf   = pred / "confidence_input_model_0.json"
    pdbs = list(pred.glob("*.pdb"))
    iptm = json.loads(sf.read_text()).get("iptm", 0.0) if sf.exists() else 0.0
    pdb  = str(pdbs[0]) if pdbs else ""
    return iptm, pdb


# ── Run ProteinMPNN ───────────────────────────────────────────────────────────

def run_mpnn():
    """Run ProteinMPNN on the FULL complex, redesigning junction residues."""
    mpnn_out = OUT_DIR / "mpnn_output"
    mpnn_out.mkdir(exist_ok=True)

    # Write fixed-positions jsonl:
    # Designable: chain A residues 79-93, chain B residues 75-89
    # Fixed: everything else
    # ProteinMPNN fixed_positions format: {"A": [1..78, 94..109], "B": [...], "C": [...], ...}
    comp_a_len = 109; comp_b_len = 143
    fixed_A = list(range(1, 79)) + list(range(94, comp_a_len + 1))
    fixed_B = list(range(1, 75)) + list(range(90, comp_b_len + 1))
    fixed_C = list(range(1, 83))
    fixed_D = list(range(1, 73))
    fixed_E = list(range(1, 126))

    fixed_jsonl = mpnn_out / "fixed_positions.jsonl"
    import json as _json
    with open(fixed_jsonl, "w") as f:
        f.write(_json.dumps({
            "complex": {
                "A": fixed_A, "B": fixed_B,
                "C": fixed_C, "D": fixed_D, "E": fixed_E
            }
        }) + "\n")

    # Copy PDB into working dir
    pdb_copy = mpnn_out / "complex.pdb"
    shutil.copy(FULL_PDB, pdb_copy)

    cmd = [
        str(MPNN_PYTHON), str(MPNN_SCRIPT),
        "--pdb_path", str(pdb_copy),
        "--out_folder", str(mpnn_out),
        "--num_seq_per_target", str(N_SEQS),
        "--sampling_temp", "0.1",
        "--fixed_positions_jsonl", str(fixed_jsonl),
        "--batch_size", "1",
    ]
    print("Running ProteinMPNN ...", flush=True)
    result = subprocess.run(cmd, capture_output=False, timeout=600)
    if result.returncode != 0:
        print("ProteinMPNN failed — check output above")
        return []

    # Parse output FASTA (seqs/complex.fa)
    fa = mpnn_out / "seqs" / "complex.fa"
    if not fa.exists():
        print(f"No FASTA at {fa}")
        return []

    entries = []
    name = seq = score = ""
    with open(fa) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    entries.append({"name": name, "seq": seq, "score": float(score) if score else 0.0})
                name = line[1:].split(",")[0]
                score_part = [p for p in line.split(",") if "score=" in p]
                score = score_part[0].split("=")[1] if score_part else "0"
                seq = ""
            else:
                seq += line
    if seq:
        entries.append({"name": name, "seq": seq, "score": float(score) if score else 0.0})

    print(f"Parsed {len(entries)} MPNN sequences")
    return entries


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    print("=== AND-gate junction optimization with ProteinMPNN ===\n")

    # Original sequences for reference
    orig_A = NLS + SOX2_BINDER + GGS3 + BRD4_FRAG_N
    orig_B = NLS + SOX9_BINDER + GGS3 + BRD4_FRAG_C
    print(f"Original CompA junction (79-93): {orig_A[78:93]}")
    print(f"Original CompB junction (75-89): {orig_B[74:89]}")

    # Run MPNN
    entries = run_mpnn()
    if not entries:
        print("No MPNN sequences produced.")
        return

    # The MPNN output FASTA encodes all 5 chains concatenated.
    # Split by chain lengths: A=109, B=143, C=82, D=72, E=125
    chain_lens = [109, 143, 82, 72, 125]
    results = []

    # Sort by MPNN score (lower = better in ProteinMPNN)
    entries.sort(key=lambda x: x["score"])

    print(f"\nTop {BOLTZ_TOP_N} by MPNN score → Boltz-2 validation:")
    for i, e in enumerate(entries[:BOLTZ_TOP_N]):
        full_seq = e["seq"]
        # Split concatenated sequence into chains
        pos = 0
        chain_seqs = []
        for l in chain_lens:
            chain_seqs.append(full_seq[pos:pos+l])
            pos += l
        comp_a, comp_b = chain_seqs[0], chain_seqs[1]

        junc_a = comp_a[78:93]
        junc_b = comp_b[74:89]
        print(f"\n[{i+1}/{BOLTZ_TOP_N}] MPNN score={e['score']:.3f}")
        print(f"  CompA junction: {junc_a}  (orig: {orig_A[78:93]})")
        print(f"  CompB junction: {junc_b}  (orig: {orig_B[74:89]})")

        wd = OUT_DIR / f"boltz_seq{i+1:02d}"
        iptm, pdb = run_boltz_5chain(comp_a, comp_b, wd)
        cf_a = brd4_cf(pdb, 'A', 'E') if pdb else 0.0
        cf_b = brd4_cf(pdb, 'B', 'E') if pdb else 0.0
        cf = max(cf_a, cf_b)
        passed = iptm >= IPTM_THRESH and cf >= CF_THRESH
        print(f"  Boltz-2: ipTM={iptm:.3f}  BRD4-cf={cf:.3f}  {'PASS' if passed else 'fail'}")

        results.append({
            "rank": i+1, "mpnn_score": e["score"],
            "comp_a": comp_a, "comp_b": comp_b,
            "junc_a": junc_a, "junc_b": junc_b,
            "iptm": iptm, "cf_brd4": cf, "passed": passed,
        })

    results.sort(key=lambda r: -(r["iptm"] + 0.3 * r["cf_brd4"]))

    print("\n\n=== Junction optimization summary ===")
    print(f"{'Rk':>3} {'ipTM':>7} {'CF':>6}  CompA junction     CompB junction")
    print("-" * 65)
    for i, r in enumerate(results, 1):
        flag = " ✓" if r["passed"] else "  "
        print(f"{i:>3} {r['iptm']:>7.3f} {r['cf_brd4']:>6.3f}{flag}  {r['junc_a']}  {r['junc_b']}")

    # Save
    tsv = PROJECT / "results" / "junction_optimization.tsv"
    with open(tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(results)
    print(f"\nSaved → {tsv}")

    # Write best FASTA
    best = results[0]
    fasta = PROJECT / "results" / "andgate_optimized.fasta"
    with open(fasta, "w") as f:
        f.write(f">ANDGATE_CompA_v1_mpnn_opt\n{best['comp_a']}\n")
        f.write(f">ANDGATE_CompB_v1_mpnn_opt\n{best['comp_b']}\n")
    print(f"Best optimized sequences → {fasta}")
    print(f"\nLead: ipTM={best['iptm']:.3f}  BRD4-cf={best['cf_brd4']:.3f}")
    print(f"  CompA: {best['comp_a']}")
    print(f"  CompB: {best['comp_b']}")


if __name__ == "__main__":
    main()

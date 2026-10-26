#!/usr/bin/env python3
"""
Boltz-2 AND-gate validation.

For each split version, run three conditions:
  FULL:  CompA + CompB + SOX2_HMG + SOX9_HMG + BRD4_BD1  (5-chain, both present)
  NO_B:  CompA          + SOX2_HMG + SOX9_HMG + BRD4_BD1  (4-chain, missing CompB)
  NO_A:          CompB  + SOX2_HMG + SOX9_HMG + BRD4_BD1  (4-chain, missing CompA)

Key metric: BRD4 contact fraction (CF) of fused BRD4-fragment residues against BRD4_BD1.
AND-gate criterion: FULL CF > 0.10  and  NO_A CF < 0.05  and  NO_B CF < 0.05
"""

import csv, json, subprocess
from pathlib import Path
import numpy as np

PROJECT = Path("/mnt/d/gbm-ANDGATE")
BOLTZ   = Path("/home/ayushregmibagale/miniconda3/envs/boltz2/bin/boltz")
OUT     = PROJECT / "validation" / "boltz_outputs"
OUT.mkdir(parents=True, exist_ok=True)

SOX2_HMG = "DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLM"
SOX9_HMG = "PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"
BRD4_BD1 = "NPPPPETSNPNKPKRQTNQLQYLLRVVLKTLWKHQFAWPFQQPVDAVKLNLPDYYKIIKTPMDMGTIKKRLENNYYWNAQECIQDFNTMFTNCYIYNKPGDDIVLMAEALEKLFLQKINELPTEE"
CONTACT_CUTOFF = 8.0


def load_constructs():
    tsv = PROJECT / "results" / "andgate_constructs.tsv"
    return {r["split"]: r for r in csv.DictReader(open(tsv), delimiter="\t")}


def write_yaml(yaml_path, chains):
    blocks = [f"  - protein:\n      id: {cid}\n      sequence: {seq}\n      msa: empty"
              for cid, seq in chains.items()]
    yaml_path.write_text("sequences:\n" + "\n".join(blocks) + "\n")


def run_boltz(work_dir, yaml_path):
    pred_dir = work_dir / "boltz_results_input" / "predictions" / "input"
    if pred_dir.exists():
        return
    subprocess.run(
        [str(BOLTZ), "predict", str(yaml_path),
         "--out_dir", str(work_dir),
         "--recycling_steps", "3", "--sampling_steps", "200",
         "--output_format", "pdb", "--no_kernels"],
        capture_output=True, timeout=1200,
    )


def parse_scores(work_dir):
    sf  = work_dir / "boltz_results_input" / "predictions" / "input" / "confidence_input_model_0.json"
    pdb = list((work_dir / "boltz_results_input" / "predictions" / "input").glob("*.pdb"))
    iptm = ptm = 0.0
    if sf.exists():
        d = json.load(open(sf))
        iptm, ptm = d.get("iptm", 0.0), d.get("ptm", 0.0)
    return iptm, ptm, str(pdb[0]) if pdb else ""


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


def brd4_cf(pdb_path, frag_chain, brd4_chain):
    ca_f = parse_ca(pdb_path, frag_chain)
    ca_b = parse_ca(pdb_path, brd4_chain)
    if not ca_f or not ca_b:
        return 0.0
    brd4_arr = np.array(list(ca_b.values()))
    return sum(np.linalg.norm(brd4_arr - xyz, axis=1).min() < CONTACT_CUTOFF
               for xyz in ca_f.values()) / len(ca_f)


def run_split(split_name, row):
    comp_a, comp_b = row["comp_a"], row["comp_b"]
    configs = {
        "FULL": {"A": comp_a, "B": comp_b, "C": SOX2_HMG, "D": SOX9_HMG, "E": BRD4_BD1},
        "NO_B": {"A": comp_a,              "C": SOX2_HMG, "D": SOX9_HMG, "E": BRD4_BD1},
        "NO_A": {              "B": comp_b, "C": SOX2_HMG, "D": SOX9_HMG, "E": BRD4_BD1},
    }
    # Boltz renumbers chains sequentially in the output PDB
    chain_maps = {
        "FULL": {"CompA":"A","CompB":"B","BRD4":"E"},
        "NO_B": {"CompA":"A","BRD4":"D"},
        "NO_A": {"CompB":"A","BRD4":"D"},
    }
    results = {}
    for cond, chains in configs.items():
        wd = OUT / f"{split_name}_{cond}"
        yp = wd / "input.yaml"
        wd.mkdir(parents=True, exist_ok=True)
        write_yaml(yp, chains)
        print(f"  Boltz-2 {split_name} {cond} ...", flush=True)
        run_boltz(wd, yp)
        iptm, ptm, pdb = parse_scores(wd)
        cmap = chain_maps[cond]
        cf_a = brd4_cf(pdb, cmap["CompA"], cmap["BRD4"]) if "CompA" in cmap else 0.0
        cf_b = brd4_cf(pdb, cmap["CompB"], cmap["BRD4"]) if "CompB" in cmap else 0.0
        cf = max(cf_a, cf_b)
        results[cond] = {"iptm": iptm, "cf_brd4": cf, "pdb": pdb}
        print(f"    → ipTM={iptm:.3f}  BRD4-cf={cf:.3f}")
    return results


def main():
    constructs = load_constructs()
    all_rows = []
    for sp in ["v1", "v2", "v3"]:
        print(f"\n── {sp} ──")
        res = run_split(sp, constructs[sp])
        for cond, r in res.items():
            all_rows.append({"split": sp, "cond": cond,
                             "iptm": r["iptm"], "cf_brd4": r["cf_brd4"]})

    print("\n=== Summary ===")
    print(f"{'Split':>6} {'Cond':>6} {'ipTM':>7} {'BRD4-CF':>8}")
    print("-" * 38)
    for r in all_rows:
        flag = "  ← AND-gate" if r["cond"]=="FULL" and r["cf_brd4"]>0.10 else ""
        print(f"{r['split']:>6} {r['cond']:>6} {r['iptm']:>7.3f} {r['cf_brd4']:>8.3f}{flag}")

    tsv = PROJECT / "results" / "andgate_validation.tsv"
    with open(tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["split","cond","iptm","cf_brd4"], delimiter="\t")
        w.writeheader()
        w.writerows(all_rows)
    print(f"\nSaved → {tsv}")


if __name__ == "__main__":
    main()

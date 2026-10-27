#!/usr/bin/env python3
"""
Run AND-gate MD replicates 2 and 3.

Starts from rep1 NPT equilibration, re-assigns velocities with different seeds.

Usage:
  python md/run_replicates.py           # run / resume
  python md/run_replicates.py --status  # show progress
"""

import os, sys, json, shutil, subprocess, argparse
import numpy as np
from pathlib import Path

BASE     = Path("/mnt/d/gbm-ANDGATE")
REP1_DIR = BASE / "md" / "rep1"
MDP_DIR  = BASE / "md" / "mdp"
GMX      = "/home/ayushregmibagale/software/gmx-plumed/bin/gmx"

REPLICATES = [
    {"id": "rep2", "seed": 51837},
    {"id": "rep3", "seed": 74293},
]

CONTACT_CUTOFF = 8.0
CF_THRESH      = 0.08


def get_state(rep_dir):
    s = rep_dir / "state.txt"
    return s.read_text().strip() if s.exists() else "PENDING"


def set_state(rep_dir, s):
    (rep_dir / "state.txt").write_text(s + "\n")
    print(f"  → {s}")


def gmx_run(args, log, cwd, stdin=None):
    env = os.environ.copy()
    env["GMX_MAXBACKUP"] = "-1"
    with open(log, "a") as f:
        r = subprocess.run(f"{GMX} {args}", shell=True, cwd=cwd, env=env,
                           input=stdin, text=True,
                           stdout=f, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(f"  FAILED (rc={r.returncode}) — see {log.name}")
        sys.exit(1)


def setup_rep(rep_dir, seed):
    rep_dir.mkdir(parents=True, exist_ok=True)
    for f in ["topol.top", "npt.gro", "start.pdb"]:
        shutil.copy2(REP1_DIR / f, rep_dir / f)
    for f in REP1_DIR.glob("topol_*.itp"):
        shutil.copy2(f, rep_dir / f.name)
    for f in REP1_DIR.glob("posre*.itp"):
        shutil.copy2(f, rep_dir / f.name)
    mdp = (MDP_DIR / "prod_50ns_newvel.mdp").read_text().replace("SEED_PLACEHOLDER", str(seed))
    (rep_dir / "prod_newvel.mdp").write_text(mdp)
    set_state(rep_dir, "SETUP_DONE")


def run_prod(rep_dir):
    print("  grompp...")
    gmx_run("grompp -f prod_newvel.mdp -c npt.gro -p topol.top -o prod.tpr -maxwarn 2",
            log=rep_dir / "grompp_prod.log", cwd=rep_dir)
    print("  mdrun (50 ns)...")
    gmx_run("mdrun -v -deffnm prod -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu",
            log=rep_dir / "mdrun_prod.log", cwd=rep_dir)
    set_state(rep_dir, "PROD_DONE")


def run_analysis(rep_dir, rep_id):
    try:
        import MDAnalysis as mda
        from MDAnalysis.lib.distances import capped_distance
    except ImportError:
        print("  MDAnalysis not found")
        sys.exit(1)

    u = mda.Universe(str(rep_dir / "prod.tpr"), str(rep_dir / "prod.xtc"))

    sox2_arm = u.select_atoms("protein and name CA and resnum 1:81")
    sox9_arm = u.select_atoms("protein and name CA and resnum 110:186")
    frag_n   = u.select_atoms("protein and name CA and resnum 91:109")
    frag_c   = u.select_atoms("protein and name CA and resnum 196:252")
    sox2_tgt = u.select_atoms("protein and name CA and resnum 253:334")
    sox9_tgt = u.select_atoms("protein and name CA and resnum 335:406")
    brd4_tgt = u.select_atoms("protein and name CA and resnum 407:531")

    cf_sox2_f, cf_sox9_f, cf_brd4_f = [], [], []
    for ts in u.trajectory:
        box = ts.dimensions
        def cf(q_sel, t_sel):
            q, t = q_sel.positions, t_sel.positions
            pairs, _ = capped_distance(q, t, max_cutoff=CONTACT_CUTOFF, box=box)
            return len(np.unique(pairs[:, 0])) / len(q) if len(pairs) and len(q) > 0 else 0.0
        cf_sox2_f.append(cf(sox2_arm, sox2_tgt))
        cf_sox9_f.append(cf(sox9_arm, sox9_tgt))
        brd4_q = np.vstack([frag_n.positions, frag_c.positions])
        pairs, _ = capped_distance(brd4_q, brd4_tgt.positions, max_cutoff=CONTACT_CUTOFF, box=box)
        cf_brd4_f.append(len(np.unique(pairs[:, 0])) / len(brd4_q) if len(pairs) else 0.0)

    cs2, cs9, cb4 = map(np.array, [cf_sox2_f, cf_sox9_f, cf_brd4_f])
    pct = float((cb4 > CF_THRESH).mean() * 100)
    metrics = {
        "rep": rep_id, "n_frames": int(len(cb4)),
        "time_ns": float(u.trajectory.n_frames * u.trajectory.dt / 1000),
        "cf_sox2_mean": float(cs2.mean()), "cf_sox9_mean": float(cs9.mean()),
        "cf_brd4_mean": float(cb4.mean()), "pct_brd4_pass": pct,
        "pass": bool(pct > 50.0),
    }
    (rep_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))
    np.savetxt(rep_dir / "contact_fractions.txt",
               np.column_stack([np.arange(len(cb4)) * u.trajectory.dt / 1000, cs2, cs9, cb4]),
               header="time_ns  cf_sox2  cf_sox9  cf_brd4", fmt="%.3f")
    print(f"  {rep_id}: sox2={metrics['cf_sox2_mean']:.3f}  sox9={metrics['cf_sox9_mean']:.3f}  "
          f"brd4={metrics['cf_brd4_mean']:.3f}  {pct:.1f}% pass  {'PASS' if metrics['pass'] else 'fail'}")
    set_state(rep_dir, "ANALYSIS_DONE")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--status", action="store_true")
    args = parser.parse_args()

    rep1_state = get_state(REP1_DIR)
    if rep1_state != "ANALYSIS_DONE":
        print(f"Rep1 not complete (state: {rep1_state}). Run: python md/run_andgate_md.py")
        sys.exit(1)

    for rep in REPLICATES:
        rep_dir = BASE / "md" / rep["id"]
        state   = get_state(rep_dir)
        if args.status:
            print(f"{rep['id']}: {state}")
            mj = rep_dir / "metrics.json"
            if mj.exists():
                print(json.dumps(json.loads(mj.read_text()), indent=2))
            continue

        print(f"\n=== {rep['id']} (seed={rep['seed']}) — state: {state} ===")
        if state == "PENDING":
            setup_rep(rep_dir, rep["seed"]); state = "SETUP_DONE"
        if state == "SETUP_DONE":
            print("  Production MD (50 ns)...")
            run_prod(rep_dir); state = "PROD_DONE"
        if state == "PROD_DONE":
            print("  Analysis...")
            run_analysis(rep_dir, rep["id"])

    if not args.status:
        print("\nAll replicates complete.")


if __name__ == "__main__":
    main()

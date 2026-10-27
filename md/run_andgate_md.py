#!/usr/bin/env python3
"""
GROMACS MD for AND-gate lead construct — 3 × 50 ns.

Starting structure: Boltz-2 5-chain prediction of optimized junction lead
(ProteinMPNN seq01, ipTM=0.489, BRD4-CF=0.210).

Chain layout (after unique residue renumbering):
  Chain A: CompA    (109 aa, res   1–109)  NLS–SOX2binder–GGS3–BRD4frag_N
  Chain B: CompB    (143 aa, res 110–252)  NLS–SOX9binder–GGS3–BRD4frag_C
  Chain C: SOX2 HMG ( 82 aa, res 253–334)
  Chain D: SOX9 HMG ( 72 aa, res 335–406)
  Chain E: BRD4 BD1 (125 aa, res 407–531)

Key sub-regions (for analysis):
  CompA_SOX2arm = res  1– 81  (NLS 1-10, SOX2binder 11-81)
  CompA_fragN   = res 91–109  (GGS3 82-90, BRD4frag_N 91-109)
  CompB_SOX9arm = res 110–186 (NLS 110-119, SOX9binder 120-186)
  CompB_fragC   = res 196–252 (GGS3 187-195, BRD4frag_C 196-252)

Analysis metrics:
  cf_sox2:  Cα contacts between CompA_SOX2arm and SOX2 HMG  (cutoff 8 Å)
  cf_sox9:  Cα contacts between CompB_SOX9arm and SOX9 HMG
  cf_brd4:  Cα contacts between (frag_N ∪ frag_C) and BRD4 BD1  ← AND-gate readout
Pass criterion: cf_brd4 > 0.08 in > 50% of frames (reconstituted binder stays engaged).

Pipeline (resumable via state.txt):
  PENDING → STRUCTURE_DONE → TOPOLOGY_DONE → EM_DONE
          → NVT_DONE → NPT_DONE → PROD_DONE → ANALYSIS_DONE

Usage:
  python md/run_andgate_md.py           # run / resume rep1
  python md/run_andgate_md.py --status  # show current stage
"""

import os, sys, json, subprocess, argparse
import numpy as np
from pathlib import Path

BASE    = Path("/mnt/d/gbm-ANDGATE")
MD_DIR  = BASE / "md" / "rep1"
MDP_DIR = BASE / "md" / "mdp"
SRC_PDB = (BASE / "design" / "mpnn_junctions" / "boltz_seq01"
           / "boltz_results_input" / "predictions" / "input" / "input_model_0.pdb")
GMX     = "/home/ayushregmibagale/software/gmx-plumed/bin/gmx"
STATE   = MD_DIR / "state.txt"

# Chain lengths → cumulative offsets for unique residue numbering
CHAIN_LENS = {"A": 109, "B": 143, "C": 82, "D": 72, "E": 125}
CONTACT_CUTOFF = 8.0
CF_THRESH      = 0.08


def get_state():
    return STATE.read_text().strip() if STATE.exists() else "PENDING"


def set_state(s):
    STATE.write_text(s + "\n")
    print(f"  → state: {s}")


def gmx(args, log, stdin=None):
    env = os.environ.copy()
    env["GMX_MAXBACKUP"] = "-1"
    with open(log, "a") as f:
        r = subprocess.run(f"{GMX} {args}", shell=True, cwd=MD_DIR, env=env,
                           input=stdin, text=True,
                           stdout=f, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(f"  FAILED (rc={r.returncode}) — see {log.name}")
        sys.exit(1)


# ── Stage 1: Fix PDB ──────────────────────────────────────────────────────────
AA1_TO_3 = {
    "A":"ALA","C":"CYS","D":"ASP","E":"GLU","F":"PHE","G":"GLY","H":"HIS",
    "I":"ILE","K":"LYS","L":"LEU","M":"MET","N":"ASN","P":"PRO","Q":"GLN",
    "R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR",
}

# Per-chain sequences (1-indexed) for resolving UNK residues
CHAIN_SEQS = {
    "A": "PKKKRKVGGSPYEIPEEDKEWDEYSFAMVEVHERLMEAERQLRLLSPEEEEKKRELEEEIERARKELEEVHKRGEAIERKKAGGGGGAAENSASIEELLEELRKKFKEN",
    # First residue of B is UNK with only backbone+CB in Boltz-2 output — use ALA
    "B": "AKKKRKVGGSMEEVRELLERARELEEKVRELEDDPSKREERERLLEEAERLREEALDLLVEQIREPLLQWMREEGLLVGSGPGGEPTEENKEKMLKELEERSEQYKEAALKFGESEEKAEEVKEYSLKLNEQVKEGWKANNQK",
}


def stage_structure():
    """Renumber residues to be globally unique; add TER between chains.
    Fixes UNK residue labels using known chain sequences (Boltz-2 artifact)."""
    MD_DIR.mkdir(parents=True, exist_ok=True)
    offsets = {}
    cum = 0
    for chain, length in CHAIN_LENS.items():
        offsets[chain] = cum
        cum += length

    lines = SRC_PDB.read_text().splitlines()
    fixed = []
    prev_chain = None
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            chain = line[21]
            if prev_chain is not None and chain != prev_chain:
                fixed.append("TER")
            if chain in offsets:
                orig = int(line[22:26])
                new  = orig + offsets[chain]
                line = line[:22] + f"{new:4d}" + line[26:]
                # Fix UNK residue names using known sequences
                if line[17:20].strip() == "UNK" and chain in CHAIN_SEQS:
                    aa1 = CHAIN_SEQS[chain][orig - 1]  # orig is 1-indexed
                    aa3 = AA1_TO_3.get(aa1, "UNK")
                    line = line[:17] + f"{aa3:<3s}" + line[20:]
            prev_chain = chain
        fixed.append(line)
    fixed.append("TER")
    fixed.append("END")

    (MD_DIR / "start.pdb").write_text("\n".join(fixed) + "\n")
    unk_remaining = sum(1 for l in fixed if "UNK" in l and l.startswith("ATOM"))
    print(f"  Written start.pdb — {cum} residues across 5 chains (UNK remaining: {unk_remaining})")
    set_state("STRUCTURE_DONE")


# ── Stage 2: Topology ─────────────────────────────────────────────────────────
def stage_topology():
    print("  pdb2gmx (AMBER99SB-ILDN, TIP3P)...")
    gmx("pdb2gmx -f start.pdb -o protein.gro -p topol.top -i posre.itp "
        "-ff amber99sb-ildn -water tip3p -ignh -missing",
        log=MD_DIR / "pdb2gmx.log")

    print("  editconf (dodecahedron, 1.2 nm)...")
    gmx("editconf -f protein.gro -o boxed.gro -bt dodecahedron -d 1.2",
        log=MD_DIR / "editconf.log")

    print("  solvate...")
    gmx("solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top",
        log=MD_DIR / "solvate.log")

    print("  genion (0.15 M NaCl)...")
    gmx(f"grompp -f {MDP_DIR}/em.mdp -c solvated.gro -p topol.top "
        "-o ions.tpr -maxwarn 2",
        log=MD_DIR / "grompp_ions.log")

    r = subprocess.run(f"{GMX} make_ndx -f solvated.gro -o /dev/null",
                       shell=True, cwd=MD_DIR, input="q\n", text=True,
                       capture_output=True, env={**os.environ, "GMX_MAXBACKUP": "-1"})
    sol_group = "15"
    for ln in (r.stdout + r.stderr).split("\n"):
        if "SOL" in ln and ln.strip() and ln.strip()[0].isdigit():
            sol_group = ln.strip().split()[0]; break

    gmx(f"genion -s ions.tpr -o ions.gro -p topol.top "
        "-pname NA -nname CL -neutral -conc 0.15",
        log=MD_DIR / "genion.log", stdin=f"{sol_group}\n")

    set_state("TOPOLOGY_DONE")


# ── Stage 3: EM ───────────────────────────────────────────────────────────────
def stage_em():
    print("  Energy minimisation...")
    gmx(f"grompp -f {MDP_DIR}/em.mdp -c ions.gro -p topol.top -o em.tpr -maxwarn 2",
        log=MD_DIR / "grompp_em.log")
    gmx("mdrun -v -deffnm em -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu",
        log=MD_DIR / "mdrun_em.log")
    set_state("EM_DONE")


# ── Stage 4: NVT ──────────────────────────────────────────────────────────────
def stage_nvt():
    print("  NVT equilibration (200 ps)...")
    gmx(f"grompp -f {MDP_DIR}/nvt.mdp -c em.gro -r em.gro -p topol.top "
        "-o nvt.tpr -maxwarn 2",
        log=MD_DIR / "grompp_nvt.log")
    gmx("mdrun -v -deffnm nvt -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu",
        log=MD_DIR / "mdrun_nvt.log")
    set_state("NVT_DONE")


# ── Stage 5: NPT ──────────────────────────────────────────────────────────────
def stage_npt():
    print("  NPT equilibration (500 ps)...")
    gmx(f"grompp -f {MDP_DIR}/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt "
        "-p topol.top -o npt.tpr -maxwarn 2",
        log=MD_DIR / "grompp_npt.log")
    gmx("mdrun -v -deffnm npt -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu",
        log=MD_DIR / "mdrun_npt.log")
    set_state("NPT_DONE")


# ── Stage 6: Production ───────────────────────────────────────────────────────
def stage_prod():
    print("  Production MD (50 ns)...")
    gmx(f"grompp -f {MDP_DIR}/prod_50ns.mdp -c npt.gro -t npt.cpt "
        "-p topol.top -o prod.tpr -maxwarn 2",
        log=MD_DIR / "grompp_prod.log")
    gmx("mdrun -v -deffnm prod -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu",
        log=MD_DIR / "mdrun_prod.log")
    set_state("PROD_DONE")


# ── Stage 7: Analysis ─────────────────────────────────────────────────────────
def stage_analysis():
    """Track contact fractions: SOX2-arm↔SOX2, SOX9-arm↔SOX9, (fragN+fragC)↔BRD4."""
    try:
        import MDAnalysis as mda
        from MDAnalysis.lib.distances import capped_distance
    except ImportError:
        print("  MDAnalysis not found — conda install -c conda-forge mdanalysis")
        sys.exit(1)

    print("  Loading trajectory...")
    u = mda.Universe(str(MD_DIR / "prod.tpr"), str(MD_DIR / "prod.xtc"))

    # Sub-region selections (renumbered residues)
    sox2_arm = u.select_atoms("protein and name CA and resnum 1:81")
    sox9_arm = u.select_atoms("protein and name CA and resnum 110:186")
    frag_n   = u.select_atoms("protein and name CA and resnum 91:109")
    frag_c   = u.select_atoms("protein and name CA and resnum 196:252")
    sox2_tgt = u.select_atoms("protein and name CA and resnum 253:334")
    sox9_tgt = u.select_atoms("protein and name CA and resnum 335:406")
    brd4_tgt = u.select_atoms("protein and name CA and resnum 407:531")

    print(f"  Atoms — sox2_arm:{len(sox2_arm)} sox9_arm:{len(sox9_arm)} "
          f"frag_n:{len(frag_n)} frag_c:{len(frag_c)} brd4:{len(brd4_tgt)}")

    cf_sox2_f, cf_sox9_f, cf_brd4_f = [], [], []

    for ts in u.trajectory:
        box = ts.dimensions

        def cf(query, target):
            q, t = query.positions, target.positions
            if len(q) == 0 or len(t) == 0:
                return 0.0
            pairs, _ = capped_distance(q, t, max_cutoff=CONTACT_CUTOFF, box=box)
            return len(np.unique(pairs[:, 0])) / len(q) if len(pairs) else 0.0

        cf_sox2_f.append(cf(sox2_arm, sox2_tgt))
        cf_sox9_f.append(cf(sox9_arm, sox9_tgt))

        # BRD4 contact: combined frag_N + frag_C vs BRD4
        brd4_query_pos = np.vstack([frag_n.positions, frag_c.positions])
        n_query = len(brd4_query_pos)
        t = brd4_tgt.positions
        pairs, _ = capped_distance(brd4_query_pos, t, max_cutoff=CONTACT_CUTOFF, box=box)
        brd4_cf = len(np.unique(pairs[:, 0])) / n_query if len(pairs) and n_query > 0 else 0.0
        cf_brd4_f.append(brd4_cf)

    cs2 = np.array(cf_sox2_f)
    cs9 = np.array(cf_sox9_f)
    cb4 = np.array(cf_brd4_f)
    pct_pass = float((cb4 > CF_THRESH).mean() * 100)

    metrics = {
        "n_frames":      int(len(cb4)),
        "time_ns":       float(u.trajectory.n_frames * u.trajectory.dt / 1000),
        "cf_sox2_mean":  float(cs2.mean()),
        "cf_sox9_mean":  float(cs9.mean()),
        "cf_brd4_mean":  float(cb4.mean()),
        "pct_brd4_pass": pct_pass,
        "pass":          bool(pct_pass > 50.0),
    }
    (MD_DIR / "metrics.json").write_text(json.dumps(metrics, indent=2))

    np.savetxt(MD_DIR / "contact_fractions.txt",
               np.column_stack([np.arange(len(cb4)) * u.trajectory.dt / 1000, cs2, cs9, cb4]),
               header="time_ns  cf_sox2  cf_sox9  cf_brd4", fmt="%.3f")

    print(f"\n{'='*60}")
    print(f"  AND-GATE MD RESULTS (rep1)")
    print(f"  Frames: {metrics['n_frames']}  ({metrics['time_ns']:.1f} ns)")
    print(f"  CF sox2-arm↔SOX2: {metrics['cf_sox2_mean']:.3f}")
    print(f"  CF sox9-arm↔SOX9: {metrics['cf_sox9_mean']:.3f}")
    print(f"  CF (fragN+fragC)↔BRD4: {metrics['cf_brd4_mean']:.3f}")
    print(f"  % frames BRD4 CF > 0.08: {pct_pass:.1f}%")
    print(f"  PASS (>50%): {metrics['pass']}")
    print(f"{'='*60}")

    set_state("ANALYSIS_DONE")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--status", action="store_true")
    args = parser.parse_args()

    state = get_state()
    if args.status:
        print(f"State: {state}")
        mj = MD_DIR / "metrics.json"
        if mj.exists():
            print(json.dumps(json.loads(mj.read_text()), indent=2))
        return

    print(f"AND-gate MD — rep1  (state: {state})")

    steps = [
        ("PENDING",        "[1/7] Fixing PDB...",           stage_structure),
        ("STRUCTURE_DONE", "[2/7] Building topology...",     stage_topology),
        ("TOPOLOGY_DONE",  "[3/7] Energy minimisation...",   stage_em),
        ("EM_DONE",        "[4/7] NVT equilibration...",     stage_nvt),
        ("NVT_DONE",       "[5/7] NPT equilibration...",     stage_npt),
        ("NPT_DONE",       "[6/7] Production MD (50 ns)...", stage_prod),
        ("PROD_DONE",      "[7/7] Analysis...",              stage_analysis),
    ]

    for from_state, label, fn in steps:
        if get_state() == from_state:
            print(f"\n{label}")
            fn()

    if get_state() == "ANALYSIS_DONE":
        print("\nAll stages complete. Run reps 2+3 with: python md/run_replicates.py")


if __name__ == "__main__":
    main()

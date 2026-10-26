#!/usr/bin/env python3
"""
Build AND-gate two-component constructs from binder sequences.

Produces andgate_constructs.tsv with Component A and Component B sequences
for each tested BRD4 split site.
"""

import csv
from pathlib import Path

PROJECT = Path("/mnt/d/gbm-ANDGATE")
OUT     = PROJECT / "results"
OUT.mkdir(exist_ok=True)

NLS  = "PKKKRKVGGS"
GGS3 = "GGSGGSGGS"

# Binders from gbm-RIPTAC / gbm-protein-switch
SOX2_BINDER = "PYEIPEEDKEWDEYSFAMVEVHERLMEAERQLRLLSPEEEEKKRELEEEIERARKELEEVHKRGEAIERKM"   # mpnn17, 71 aa
SOX9_BINDER = "MEEVRELLERARELEEKVRELEDDPSKREERERLLEEAERLREEALDLLVEQIREPLLQWMREEGLL"         # ph_sox9_run18, 67 aa
BRD4_FULL   = "MNKSIEELLEELRKKFKENTEENKEKMLKELEERSEQYKEAALKFGESEEKAEEVKEYSLKLNEQVKEGWKANNQK"  # ph_r49c2, 76 aa

# Split sites: (name, position, rationale)
SPLIT_SITES = [
    ("v1", 19, "natural low-pLDDT loop — winner (ΔBRD4-CF=+0.203)"),
    ("v2", 38, "balanced midpoint"),
    ("v3", 46, "after Gly turn"),
]

rows = []
for name, pos, note in SPLIT_SITES:
    frag_n = BRD4_FULL[:pos]
    frag_c = BRD4_FULL[pos:]
    comp_a = NLS + SOX2_BINDER + GGS3 + frag_n
    comp_b = NLS + SOX9_BINDER + GGS3 + frag_c
    rows.append({
        "split": name, "pos": pos, "note": note,
        "frag_n": frag_n, "frag_c": frag_c,
        "comp_a": comp_a, "comp_b": comp_b,
        "len_a": len(comp_a), "len_b": len(comp_b),
    })
    print(f"\n── {name} (split after pos {pos}) — {note}")
    print(f"  frag_N: {frag_n}  ({len(frag_n)} aa)")
    print(f"  frag_C: {frag_c}  ({len(frag_c)} aa)")
    print(f"  CompA ({len(comp_a)} aa): {comp_a}")
    print(f"  CompB ({len(comp_b)} aa): {comp_b}")

tsv = OUT / "andgate_constructs.tsv"
with open(tsv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print(f"\nSaved → {tsv}")

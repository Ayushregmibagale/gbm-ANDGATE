# gbm-ANDGATE Progress Log

## Project Summary
A two-component AND-gate system that sequesters BRD4 from super-enhancers **only when both SOX2 and SOX9 are co-present** — the GSC-specific state. Normal cells expressing only one factor are protected.

Spun out from **gbm-RIPTAC** (Route 1 AND-gate extension), 2026-10-26.

---

## Current Status (2026-10-28) — MD Complete ✅

### ✅ Computational validation complete

**Binder sequences (imported from gbm-RIPTAC / gbm-protein-switch):**
| Arm | Source | ipTM | Len | Sequence |
|-----|--------|------|-----|----------|
| SOX2 binder | mpnn17 (gbm-protein-switch) | 0.897 | 71 | `PYEIPEEDKEWDEYSFAMVEVHERLMEAERQLRLLSPEEEEKKRELEEEIERARKELEEVHKRGEAIERKM` |
| SOX9 binder | ph_sox9_run18 (gbm-protein-switch) | 0.856 | 67 | `MEEVRELLERARELEEKVRELEDDPSKREERERLLEEAERLREEALDLLVEQIREPLLQWMREEGLL` |
| BRD4 binder | ph_r49c2 (gbm-RIPTAC) | 0.948 | 76 | `MNKSIEELLEELRKKFKENTEENKEKMLKELEERSEQYKEAALKFGESEEKAEEVKEYSLKLNEQVKEGWKANNQK` |

**Split site analysis:**
- pLDDT profiling of BRD4 binder structure identified a natural low-confidence loop at residues 18–24
- Split site selected at position 19 (`...ELRKKFKEN | TEENKEK...`)
- Rationale: frag_N (19 aa) is too short to independently contact BRD4; frag_C (57 aa) lacks the anchoring helix

**Boltz-2 AND-gate validation (3 split sites × 3 conditions = 9 jobs):**

| Split | Condition | ipTM | BRD4-CF | AND-gate |
|-------|-----------|------|---------|---------|
| **v1 (pos 19)** | FULL (both) | 0.386 | **0.203** | ✓ |
| v1 | NO CompB | 0.264 | **0.000** | — |
| v1 | NO CompA | 0.549 | **0.000** | — |
| v2 (pos 38) | FULL | 0.304 | 0.125 | ✓ (leakage v2-NO_B=0.031) |
| v3 (pos 46) | FULL | 0.395 | 0.162 | ✓ (leakage v3-NO_B=0.059) |

**Winner: v1 (split at pos 19)** — clean AND-gate, zero BRD4 contact when either component absent.

### GROMACS MD — 3 × 50 ns ✅ ALL PASS

5-chain complex (CompA + CompB + SOX2 + SOX9 + BRD4), full reconstituted AND-gate.

| Rep | SOX2 CF | SOX9 CF | BRD4 CF | % frames BRD4 pass |
|-----|---------|---------|---------|---------------------|
| rep1 | 0.189 | 0.162 | 0.430 | **100.0%** |
| rep2 | 0.224 | 0.153 | 0.419 | **100.0%** |
| rep3 | 0.166 | 0.144 | 0.366 | **100.0%** |
| **Mean** | **0.193** | **0.153** | **0.405** | **100.0%** |

BRD4 reconstituted binding is stable across 100% of frames in all three replicates. **PASS.**

### Final constructs (v1)

**Component A (109 aa) — `NLS–SOX2_binder–GGS3–BRD4_frag_N`**
```
PKKKRKVGGSPYEIPEEDKEWDEYSFAMVEVHERLMEAERQLRLLSPEEEEKKRELEEEIERARKELEEVHKRGEAIERKMGGSGGSGGSMNKSIEELLEELRKKFKEN
```
Breakdown: `PKKKRKVGGS` (NLS) + `PYEIPEEDKEW...RGEAIERKM` (SOX2 binder, 71 aa) + `GGSGGSGGS` (GGS3) + `MNKSIEELLEELRKKFKEN` (BRD4 frag_N, 19 aa)

**Component B (143 aa) — `NLS–SOX9_binder–GGS3–BRD4_frag_C`**
```
PKKKRKVGGSMEEVRELLERARELEEKVRELEDDPSKREERERLLEEAERLREEALDLLVEQIREPLLQWMREEGLLGGSGGSGGSTEENKEKMLKELEERSEQYKEAALKFGESEEKAEEVKEYSLKLNEQVKEGWKANNQK
```
Breakdown: `PKKKRKVGGS` (NLS) + `MEEVRELLERARELEEKVR...WREEGLL` (SOX9 binder, 67 aa) + `GGSGGSGGS` (GGS3) + `TEENKEKMLK...ANNQK` (BRD4 frag_C, 57 aa)

---

## Mechanism

```
GSC (SOX2 + SOX9 co-present):
  CompA ──→ binds SOX2 ─┐
                         ├─ SOX2:SOX9 complex brings A+B together
  CompB ──→ binds SOX9 ─┘
                         → frag_N + frag_C reconstitute BRD4 binder
                         → BRD4 sequestered from super-enhancers
                         → GSC dies ✓

Normal cell (SOX2 only or SOX9 only):
  CompA binds SOX2, frag_N (19 aa) floats → cannot contact BRD4 ✗
  CompB absent or not co-localized         → frag_C floats alone  ✗
  BRD4 activity preserved ✓
```

---

## Aims

1. ✅ **Aim 1 — Computational validation** (done 2026-10-26)
2. ✅ **MD simulation** — 3×50 ns GROMACS, BRD4 CF=0.405 at 100% frames (done 2026-10-28)
3. ⚠️ **Aim 2 — Junction optimization** — see findings below (done 2026-10-28)
4. ⏳ **Aim 3 — Expression** — HEK293T expression, nuclear localization, SOX co-IP
5. ⏳ **Aim 4 — Reconstitution assay** — Co-expression in GSC vs NPC, BRD4 ChIP, viability
6. ⏳ **Aim 5 — Orthogonality** — BRD4 function in non-GSC contexts

---

## File Map

```
gbm-ANDGATE/
├── aims.md
├── PROGRESS.md
├── structures/
│   ├── SOX2/mpnn17_sox2_binder.pdb
│   ├── SOX9/ph_sox9_run18_binder.pdb
│   └── BRD4/2OSS_BD1.pdb + ph_r49c2_brd4_binder.pdb
├── design/
│   └── build_andgate.py          Generate component sequences (all split variants)
├── validation/
│   └── validate_andgate.py       Boltz-2 FULL/NO_A/NO_B validation
├── results/
│   ├── andgate_final.fasta        Final CompA + CompB sequences (v1)
│   ├── andgate_constructs.tsv     All split site constructs
│   └── andgate_validation.tsv    Boltz-2 results (9 jobs)
└── md/
    ├── run_andgate_md.py          Rep1 pipeline (resumable)
    ├── run_replicates.py          Rep2+3
    ├── results/                   Tracked metrics (rep1–3)
    └── mdp/                      GROMACS parameter files (50 ns)
```

---

## Next Steps

### Aim 2 — ProteinMPNN Junction Optimization (2026-10-28)

Ran ProteinMPNN on the 5-chain FULL complex (20 junction designs), then validated all 20 with
Boltz-2 FULL + NO_B + NO_A conditions.

**Result: all top-5 designs by FULL BRD4-CF fail AND-gate validation.**

| Design | FULL CF | NO_B CF | NO_A CF | Result |
|--------|---------|---------|---------|--------|
| seq08 | 0.487 | 0.421 | 0.439 | fail |
| seq12 | 0.461 | 0.474 | 0.474 | fail |
| seq14 | 0.461 | 0.316 | 0.386 | fail |
| seq20 | 0.461 | 0.579 | 0.439 | fail |
| seq01 | 0.447 | 0.421 | 0.456 | fail |

**Root cause**: the MPNN fixed positions included BRD4-fragment N-terminal residues (frag_N pos
91-93 of CompA; frag_C pos 87-89 of CompB). Redesigning these residues improved reconstituted
complex assembly but inadvertently stabilized individual fragment BRD4 binding → leakage.
The original v1 split design has 0.000 leakage because the natural fragment N-termini are
unfavourable for BRD4 contact in isolation.

**Fix attempted (Aim 2b)**: constrained MPNN — fix all BRD4-fragment residues (A91-109 and B87-143);
only redesign the GGS3 linker and binder C-terminus (A78-90, B75-86). See Aim 2b below.

**Decision**: proceed to Phase 1 lab with **original v1 design** (zero leakage,
BRD4-CF=0.405 at 100% of MD frames). Aim 2b constrained MPNN attempted and also failed — see below.

### Aim 2b — Constrained MPNN (GGS3 linker only) (2026-10-28)

Re-ran ProteinMPNN with BRD4-fragment residues fully fixed (A91-109, B87-143); only the
GGS3 linker + binder C-terminus redesigned (A78-90, B75-86 — 13 and 12 residues respectively).

**Result: all top-5 designs still fail AND-gate validation.**

| Design | FULL CF | NO_B CF | NO_A CF | Result |
|--------|---------|---------|---------|--------|
| seq20 | 0.487 | 0.368 | 0.404 | fail |
| seq09 | 0.461 | 0.474 | 0.509 | fail |
| seq10 | 0.461 | 0.474 | 0.474 | fail |
| seq12 | 0.461 | 0.526 | 0.456 | fail |
| seq14 | 0.461 | 0.474 | 0.386 | fail |

**Root cause confirmed**: MPNN replaces the flexible GGS3 linker (GGSGGSGGS) with structured
variants (e.g. seq20 linker: `SGGGGGSE`). The original GGSGGSGGS composition keeps frag_N/frag_C
disordered in isolation — critical for AND-gate specificity. Any MPNN redesign of the GGS3
region breaks this property.

**Final conclusion**: The GGSGGSGGS linker sequence is load-bearing for AND-gate logic and
**must not be redesigned**. Junction optimization via MPNN is incompatible with the AND-gate
design as implemented. Aim 2 is closed; computational optimization is complete.

**Lead design**: original v1 (zero leakage, BRD4-CF=0.405 at 100% of MD frames). Ready for lab.

---

## Next Steps

1. ✅ **MD simulation**: 3×50 ns complete, BRD4 CF=0.405 at 100% frames
2. ✅ **Aim 2 / Aim 2b**: MPNN junction optimization exhausted — original v1 is optimal
3. **Codon optimization** for mammalian expression (IDT codon tool or python-codon-tables)
4. **Order synthetic genes** — two separate expression constructs (separate plasmids, same promoter)
5. **Phase 1 lab validation** (Aim 3)

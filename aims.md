# gbm-ANDGATE — Project Aims

## The problem with always-on BRD4 sequestration

The RIPTAC approach (fused SOX2 + SOX9 + BRD4 binder) sequesters BRD4 from super-enhancers in
glioblastoma stem cells (GSCs). It works because SOX2 and SOX9 are co-expressed in GSCs, and their
complex drives the oncogenic transcriptional program together with BRD4.

**The selectivity gap:** A single fused protein has all three arms always active. A cell expressing
only SOX2 (e.g., normal neural stem cells) or only SOX9 (e.g., chondrocytes) would have its BRD4
sequestered unnecessarily. This is an off-target risk.

## The AND-gate solution

Split the BRD4 binder into two non-functional halves. Distribute them across two separate proteins,
each carrying one SOX-family binder:

```
Component A:  NLS – [SOX2 binder] – [BRD4 frag N]
Component B:  NLS – [SOX9 binder] – [BRD4 frag C]
```

**Mechanism:**
- In cells with SOX2 only: Component A binds SOX2, but frag_N (19 aa) alone cannot contact BRD4.
- In cells with SOX9 only: Component B binds SOX9, but frag_C (57 aa) alone cannot contact BRD4.
- In GSCs (SOX2 + SOX9 co-expressed): SOX2 and SOX9 form their natural heteromeric complex on
  super-enhancers. This brings Component A and Component B into proximity. The two BRD4 halves
  reconstitute the full 76 aa binder → BRD4 is sequestered. GSCs lose super-enhancer activity and die.

The SOX2:SOX9 protein–protein interaction is itself the logic gate.

## Aims

1. **Aim 1 — Computational validation**: Confirm AND-gate logic with Boltz-2 structure prediction
   (BRD4 contact fraction: FULL >> NO_A == NO_B ≈ 0).

2. **Aim 2 — Junction optimization**: Use ProteinMPNN to redesign the GGS3–frag junction sequences
   in both components for improved folding and reconstitution efficiency.

3. **Aim 3 — Wet-lab expression**: Express Component A and Component B individually in HEK293T.
   Confirm nuclear localization (anti-NLS staining), SOX2/SOX9 binding (co-IP), and lack of
   independent BRD4 binding (BRD4 co-IP with single components in SOX2- or SOX9-only cells).

4. **Aim 4 — Reconstitution assay**: Co-express both components in GSC lines (SOX2+SOX9+) vs
   normal neural progenitors (SOX2-only). Measure BRD4 ChIP at super-enhancers and GSC viability.

5. **Aim 5 — Orthogonality**: Ensure neither component disrupts BRD4 function in non-GSC contexts.

## Relationship to gbm-RIPTAC

This project is a **conditional-logic extension** of gbm-RIPTAC (Route 1). The binder sequences
are taken directly from gbm-RIPTAC:
- SOX2 binder: `mpnn17` (ipTM=0.897, from gbm-protein-switch)
- SOX9 binder: `ph_sox9_run18` (ipTM=0.856, from gbm-protein-switch)
- BRD4 binder: `ph_r49c2` (ipTM=0.948, Protein Hunter hallucination in gbm-RIPTAC)

The split site (after residue 19 of BRD4 binder, natural low-pLDDT loop) was identified by
pLDDT analysis and confirmed by Boltz-2 AND-gate validation (ΔBRD4-CF = +0.203).

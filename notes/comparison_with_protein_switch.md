# gbm-ANDGATE vs gbm-protein-switch: Project Comparison

*Written 2026-10-28*

## Relationship

gbm-ANDGATE was spun off from gbm-protein-switch Track 1. It re-uses the exact same SOX2
and SOX9 binders and tests whether a non-apoptotic actuator (BRD4 sequestration) can achieve
the same GSC selectivity with a potentially cleaner safety profile.

---

## Shared core — identical in both projects

| Component | gbm-protein-switch | gbm-ANDGATE |
|-----------|-------------------|-------------|
| SOX2 binder | mpnn17 (71 aa) | mpnn17 (71 aa) |
| SOX9 binder | ph_sox9_run18 (67 aa) | ph_sox9_run18 (67 aa) |
| Linker | GGS3 (GGSGGSGGS) | GGS3 (GGSGGSGGS) |
| NLS | SV40 PKKKRKVGGS | SV40 PKKKRKVGGS |
| Architecture | NLS–SOX2_binder–GGS3–[half-A] + NLS–SOX9_binder–GGS3–[half-B] | identical |
| Sensing logic | SOX2 AND SOX9 co-present → halves reconstitute | identical |

The entire sensing apparatus is the same. Both projects are the **same two-component system
with a different payload**.

---

## What differs: the actuator

| | gbm-protein-switch | gbm-ANDGATE |
|--|--|--|
| **Payload** | Split caspase-9 (CARD-deleted, D315/G316) | Split BRD4-binder (ph_r49c2, split pos 19) |
| **Mechanism of kill** | Reconstituted caspase-9 → caspase-3/7 → apoptosis | Reconstituted BRD4-binder → sequesters BRD4 from super-enhancers → transcriptional collapse |
| **Speed** | Hours (direct apoptosis) | Days (transcriptional) |
| **Signal amplification** | Yes — one reconstituted dimer activates thousands of downstream caspases | No — stoichiometric BRD4 occupancy required |
| **Fragment sizes** | N-half 176 aa + C-half 101 aa | frag_N 19 aa + frag_C 57 aa |
| **MD metric** | Inter-half C-term distance < 15 Å | BRD4 contact fraction > 0.08 |
| **MD result** | 81.4% frames PASS (3 × 20 ns) | 100% frames PASS (3 × 50 ns) |

---

## Computational status (as of 2026-10-28)

### gbm-protein-switch ✅
- Track 1 lead: mpnn17 × ph_sox9_run18 — 3 × 20 ns MD, 81.4% frames inter-half < 15 Å
- Track 2 lead: bridge_48_s01 — single bridging miniprotein, 99.4% dual-contact

### gbm-ANDGATE ✅
- Lead: original v1 split (pos 19) — 3 × 50 ns MD, BRD4 CF=0.405 at 100% of frames
- Boltz-2 AND-gate: zero leakage in NO_A / NO_B single-component controls
- Aim 2 / 2b MPNN junction optimization: exhausted — GGS3 linker is load-bearing; original v1 optimal

Both computationally complete. Both ready for lab.

---

## Lab strategy: shared experiments

Because the sensing apparatus is identical, several experiments validate both projects at once:

| Experiment | Applies to |
|-----------|-----------|
| His-tag expression of CompA + CompB (HEK293T) | both |
| NLS staining — confirm nuclear localization | both |
| SOX2 co-IP with CompA | both |
| SOX9 co-IP with CompB | both |
| GSC vs NPC co-expression — confirm co-localization | both |

Projects diverge only at the payload assay:
- **gbm-protein-switch**: caspase-3/7 activity assay (Caspase-Glo), viability
- **gbm-ANDGATE**: BRD4 ChIP-seq (super-enhancer displacement), viability

---

## Which to prioritise?

| Criterion | gbm-protein-switch | gbm-ANDGATE |
|-----------|-------------------|-------------|
| Kill speed | Faster (hours) | Slower (days) |
| Signal amplification | Yes (catalytic) | No (stoichiometric) |
| Fragment size (cargo efficiency) | Larger (277 aa total) | Smaller (76 aa total) |
| Safety concern | Caspase-9 leakage = apoptosis | BRD4 leakage = partial transcriptional effect |
| MD confidence | 81.4% frames | 100% frames |
| AND-gate leakage test | Not yet done (MD of unbound state done: termini >25 Å, 95% frames) | Done: 0.000 leakage (Boltz-2) |

Both are strong candidates. Recommended approach: run shared sensing-apparatus experiments
in parallel, then diverge at the payload assay to compare both actuators head-to-head in GSCs.

#!/usr/bin/env bash
# Resumes ANDGATE rep3 from checkpoint after RIPTAC rep3 finishes.
GMX=/home/ayushregmibagale/software/gmx-plumed/bin/gmx
REP3=/mnt/d/gbm-ANDGATE/md/rep3
RIPTAC_STATE=/mnt/d/gbm-RIPTAC/md/rep3/state.txt

echo "[resume] Waiting for RIPTAC rep3 to finish..."
until grep -q "ANALYSIS_DONE" "$RIPTAC_STATE" 2>/dev/null; do sleep 60; done
echo "[resume] RIPTAC rep3 done. Resuming ANDGATE rep3 from checkpoint at $(date)"

export GMX_MAXBACKUP=-1
cd "$REP3"
$GMX mdrun -v -deffnm prod -cpi prod.cpt -ntmpi 1 -ntomp 8 -gpu_id 0 -nb gpu -pme gpu -bonded gpu \
    >> "$REP3/mdrun_prod.log" 2>&1

echo "[resume] ANDGATE rep3 mdrun finished at $(date)"

# Run analysis via the replicates script (rep2 already done, rep3 just finished)
cd /mnt/d/gbm-ANDGATE
python md/run_replicates.py >> md/replicates.log 2>&1
echo "[resume] ANDGATE rep3 analysis complete at $(date)"

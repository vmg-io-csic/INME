#!/bin/bash
# Resolve the repository root from the location of this script, so the script
# works from any working directory and on any machine.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Python interpreter: override with INME_PYTHON if the networks live in a venv.
PY="${INME_PYTHON:-python3}"
# Where stimuli and results are written: override with INME_OUT.
OUTROOT="${INME_OUT:-$REPO/output}"
S=/tmp/claude-1000/-home-mberraquero/c2798504-67d8-481a-a009-5fc07ff6d143/scratchpad
until grep -q "SPEEDFIX GENERADO" $S/speedfix_gen.log 2>/dev/null; do sleep 30; done
while pgrep -f run_reps.sh >/dev/null; do sleep 60; done            # esperar a la cola de repeticiones
sed -e 's#stim_mask_sweep#stim_speedfix#' -e 's#results_mask_sweep#results_speedfix#' $S/run_mask_dnn.sh > $S/run_speedfix_inner.sh
chmod +x $S/run_speedfix_inner.sh; $S/run_speedfix_inner.sh

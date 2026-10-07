#!/bin/bash
# Resolve the repository root from the location of this script, so the script
# works from any working directory and on any machine.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Python interpreter: override with INME_PYTHON if the networks live in a venv.
PY="${INME_PYTHON:-python3}"
# Where stimuli and results are written: override with INME_OUT.
OUTROOT="${INME_OUT:-$REPO/output}"
S=/tmp/claude-1000/-home-mberraquero/c2798504-67d8-481a-a009-5fc07ff6d143/scratchpad
while pgrep -f run_contrast20.sh >/dev/null; do sleep 30; done          # esperar a que termine contraste
echo "### GENERANDO"
matlab -batch "ps=parallel.Settings; ps.Pool.AutoCreate=false; run('$S/gen_reps.m')" 2>&1 | grep -E "generado|TERMINADA|Error"

cd "$REPO/dnn_inference"
for EXP in ReversePhi DH TimeReversal SecondOrder Masking; do
  ST=$OUTROOT/stimuli/$EXP; OUT=$OUTROOT/results/$EXP
  for M in raft pwcnet liteflownet2; do
    G=pwc_lfn2; [ "$M" = raft ] && G=raft
    $PY run_dnn_fullvideo_gauss.py --model $M --ckpt sintel --stimuli_dir $ST/$G/v25 --output_dir $OUT --save_flow none --gpu 0 > $S/reps_${EXP}_${M}.log 2>&1 || echo "### FALLO $EXP $M"
    echo "### FIN $EXP $M"
  done
  $PY run_multiframe_pipeline.py --model me_attention --model_dir me_attention_dual/ --checkpoint me_attention_dual/Model_example.pth.tar --stimuli_dir $ST/meattn_dual/v25 --output_dir $OUT --n_input_frames 11 --roi_half 100 --sigma_c 57.735 --tb_a 9 --tb_b 1 --skip_threshold 0 --gpu 0 > $S/reps_${EXP}_meatt.log 2>&1 || echo "### FALLO $EXP me_attention"
  echo "### FIN $EXP me_attention"
  $PY run_multiframe_pipeline.py --model dual --model_dir dual_model/ --checkpoint me_attention_dual/final_sintel_kitti.pth --stimuli_dir $ST/meattn_dual/v25 --output_dir $OUT --n_input_frames 15 --roi_half 100 --sigma_c 57.735 --tb_a 9 --tb_b 1 --skip_threshold 0 --gpu 0 > $S/reps_${EXP}_dual.log 2>&1 || echo "### FALLO $EXP dual"
  echo "### FIN $EXP dual"
  echo "### EXPERIMENT DONE: $EXP"
done
echo "### ALL DONE"

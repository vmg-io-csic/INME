#!/bin/bash
# Resolve the repository root from the location of this script, so the script
# works from any working directory and on any machine.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Python interpreter: override with INME_PYTHON if the networks live in a venv.
PY="${INME_PYTHON:-python3}"
# Where stimuli and results are written: override with INME_OUT.
OUTROOT="${INME_OUT:-$REPO/output}"

S=$OUTROOT/stimuli/mask_sweep
OUT=$OUTROOT/results/mask_sweep
cd "$REPO/dnn_inference"
for M in raft pwcnet liteflownet2; do
  G=pwc_lfn2; [ "$M" = raft ] && G=raft
  echo "### INICIO $M"
  $PY run_dnn_fullvideo_gauss.py --model $M --ckpt sintel --stimuli_dir $S/$G/v25 --output_dir $OUT --save_flow none --gpu 0 || echo "### FALLO $M"
  echo "### FIN $M"
done
echo "### INICIO me_attention"
$PY run_multiframe_pipeline.py --model me_attention --model_dir me_attention_dual/ --checkpoint me_attention_dual/Model_example.pth.tar --stimuli_dir $S/meattn_dual/v25 --output_dir $OUT --n_input_frames 11 --roi_half 100 --sigma_c 57.735 --tb_a 9 --tb_b 1 --skip_threshold 0 --gpu 0 || echo "### FALLO me_attention"
echo "### FIN me_attention"
echo "### INICIO dual"
$PY run_multiframe_pipeline.py --model dual --model_dir dual_model/ --checkpoint me_attention_dual/final_sintel_kitti.pth --stimuli_dir $S/meattn_dual/v25 --output_dir $OUT --n_input_frames 15 --roi_half 100 --sigma_c 57.735 --tb_a 9 --tb_b 1 --skip_threshold 0 --gpu 0 || echo "### FALLO dual"
echo "### FIN dual"
echo "### ALL DONE"

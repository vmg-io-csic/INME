#!/bin/bash
# Run one set of stimuli through the five networks. Usage:
#   ./run_dnn_experiment.sh <stimuli_dir> <results_dir>
# The stimuli directory holds one subdirectory per network group
# (raft, pwc_lfn2, meattn_dual), each containing a v25 folder.
set -u
# Resolve the repository root from the location of this script, so the script
# works from any working directory and on any machine.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Python interpreter: override with INME_PYTHON if the networks live in a venv.
PY="${INME_PYTHON:-python3}"
# Where stimuli and results are written: override with INME_OUT.
OUTROOT="${INME_OUT:-$REPO/output}"
ST=${1:?stimuli directory}; OUT=${2:?results directory}

cd "$REPO/dnn_inference"
for M in raft pwcnet liteflownet2; do
  G=pwc_lfn2; [ "$M" = raft ] && G=raft
  $PY run_dnn_fullvideo_gauss.py --model $M --ckpt sintel \
      --stimuli_dir "$ST/$G/v25" --output_dir "$OUT" --save_flow none --gpu 0
done
$PY run_multiframe_pipeline.py --model me_attention --model_dir me_attention_dual/ \
    --checkpoint me_attention_dual/Model_example.pth.tar --stimuli_dir "$ST/meattn_dual/v25" \
    --output_dir "$OUT" --n_input_frames 11 --roi_half 100 --sigma_c 57.735 \
    --tb_a 9 --tb_b 1 --skip_threshold 0 --gpu 0
$PY run_multiframe_pipeline.py --model dual --model_dir dual_model/ \
    --checkpoint me_attention_dual/final_sintel_kitti.pth --stimuli_dir "$ST/meattn_dual/v25" \
    --output_dir "$OUT" --n_input_frames 15 --roi_half 100 --sigma_c 57.735 \
    --tb_a 9 --tb_b 1 --skip_threshold 0 --gpu 0

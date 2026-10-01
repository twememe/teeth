#!/usr/bin/env bash
set -euo pipefail
cd /home/xinlab/home/ws/phase2
mkdir -p logs results/12_complex_prediction/glyco_sensitivity
export CUDA_VISIBLE_DEVICES=0

for route in native_blind immunogen_blind; do
  target="${route}_msa_n74_nag"
  for seed in 1 2 3; do
    out="results/12_complex_prediction/glyco_sensitivity/${route}_seed${seed}"
    result="${out}/boltz_results_${target}/predictions/${target}/confidence_${target}_model_0.json"
    if [[ -f "$result" ]]; then
      echo "SKIP existing ${route} seed ${seed}"
      continue
    fi
    echo "START $(date --iso-8601=seconds) ${route} seed ${seed}"
    mkdir -p "$out"
    venv/bin/boltz predict \
      "results/12_complex_prediction/boltz_inputs/${target}.yaml" \
      --model boltz1 --cache models/boltz --out_dir "$out" \
      --checkpoint models/boltz/boltz1_conf.ckpt --seed "$seed" \
      --accelerator gpu --devices 1 --diffusion_samples 1 \
      --recycling_steps 3 --sampling_steps 200 --no_kernels \
      --num_workers 0 --output_format pdb --write_full_pae \
      > "logs/boltz_glyco_${route}_seed${seed}.log" 2>&1
    echo "DONE  $(date --iso-8601=seconds) ${route} seed ${seed}"
  done
done
echo "GLYCO_SENSITIVITY_DONE $(date --iso-8601=seconds)"

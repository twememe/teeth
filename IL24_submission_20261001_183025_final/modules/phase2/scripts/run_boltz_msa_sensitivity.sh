#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs results/12_complex_prediction/msa_sensitivity

for route in native_blind immunogen_blind; do
  target="${route}_msa"
  for seed in 1 2 3; do
    out="results/12_complex_prediction/msa_sensitivity/${route}_seed${seed}"
    result="${out}/boltz_results_${target}/predictions/${target}/confidence_${target}_model_0.json"
    if [[ -f "$result" ]]; then
      echo "SKIP existing ${route} seed ${seed}"
      continue
    fi
    echo "RUN ${route} seed ${seed}"
    mkdir -p "$out"
    "${BOLTZ_BIN:-boltz}" predict \
      "results/12_complex_prediction/boltz_inputs/${target}.yaml" \
      --model boltz1 \
      --cache "${BOLTZ_CACHE:?Set BOLTZ_CACHE to the downloaded Boltz assets directory}" \
      --out_dir "$out" \
      --checkpoint "${BOLTZ_CACHE:?}/boltz1_conf.ckpt" \
      --seed "$seed" \
      --accelerator gpu \
      --devices 1 \
      --diffusion_samples 1 \
      --recycling_steps 3 \
      --sampling_steps 50 \
      --no_kernels \
      --num_workers 0 \
      --output_format pdb \
      --write_full_pae \
      > "logs/boltz_msa_${route}_seed${seed}.log" 2>&1
    echo "DONE ${route} seed ${seed}"
  done
done
echo MSA_SENSITIVITY_DONE

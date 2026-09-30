#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
for route in native_blind immunogen_blind; do
  for seed in $(seq 1 15); do
    out="results/12_complex_prediction/formal/${route}_seed${seed}"
    result="${out}/boltz_results_${route}/predictions/${route}/confidence_${route}_model_0.json"
    if [[ -f "$result" ]]; then continue; fi
    mkdir -p "$out"
    boltz predict "results/12_complex_prediction/boltz_inputs/${route}.yaml" \
      --model boltz1 --cache models/boltz --out_dir "$out" \
      --checkpoint models/boltz/boltz1_conf.ckpt --seed "$seed" \
      --accelerator gpu --devices 1 --diffusion_samples 1 \
      --recycling_steps 3 --sampling_steps 200 --no_kernels \
      --num_workers 0 --output_format pdb --write_full_pae \
      > "logs/boltz_formal_${route}_seed${seed}.log" 2>&1
  done
done
echo FORMAL_DONE

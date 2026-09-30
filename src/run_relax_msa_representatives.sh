#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
OUTPUT_ROOT="${PHASE2_OUTPUT_ROOT:?PHASE2_OUTPUT_ROOT is required}"
export OPENMM_CPU_THREADS=32
export OPENMM_PLATFORM=CPU
export RELAX_MAX_ITERATIONS=100
mkdir -p "$OUTPUT_ROOT/14_quality/relaxed"

for name in \
  native_blind_msa_consensus_representative_seed9 \
  immunogen_blind_msa_consensus_representative_seed7; do
  source="$OUTPUT_ROOT/14_quality/representatives_msa/${name}.pdb"
  destination="$OUTPUT_ROOT/14_quality/relaxed/${name}_relaxed.pdb"
  report="$OUTPUT_ROOT/14_quality/relaxed/${name}_relaxation.json"
  if [[ -f "$destination" ]]; then echo "SKIP $name"; continue; fi
  echo "RELAX $name"
  python3 src/relax_representative_structures.py \
    "$source" "$destination" --report "$report"
done
echo RELAX_MSA_REPRESENTATIVES_DONE

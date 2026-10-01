# Phase 2 reproduction notes

Server root: the `phase2/` directory of this repository. No absolute path is required; `bash run.sh` and `predict.py` derive it from their own location.

Environment: Python 3.10, torch 2.6.0+cu124, CUDA 12.4, NVIDIA RTX A6000, Boltz 2.2.1, PRODIGY 2.4.0.

Core commands:

```bash
cd <repository>/phase2
venv/bin/python scripts/run_imgt_numbering.py
PATH=$PWD/hmmer_bin/bin:$PATH venv/bin/python scripts/step11_antibody_structure.py
venv/bin/python scripts/generate_boltz_inputs.py
bash scripts/run_boltz_pilots.sh
bash scripts/run_boltz_formal.sh
venv/bin/python scripts/analyze_boltz_complexes.py --set formal
venv/bin/python scripts/run_prodigy_formal.py
venv/bin/python scripts/rank_formal_models.py
venv/bin/python scripts/summarize_phase2.py
```

MSA sensitivity inputs are generated locally because the server TLS path inserts an untrusted self-signed certificate for `api.colabfold.com`:

```powershell
python scripts/generate_colabfold_msa.py
python scripts/generate_boltz_inputs.py
```

After uploading `data/msa/` and the generated YAML files to the same relative paths on the server:

```bash
bash scripts/run_boltz_msa_sensitivity.sh
venv/bin/python scripts/analyze_boltz_complexes.py --set msa_sensitivity
venv/bin/python scripts/summarize_msa_sensitivity.py
```

Formal MSA, single-NAG and relaxation extensions:

```bash
bash scripts/run_boltz_formal_msa.sh
venv/bin/python scripts/analyze_boltz_complexes.py --set formal_msa
venv/bin/python scripts/run_prodigy_formal.py --set formal_msa
venv/bin/python scripts/summarize_formal_msa.py
bash scripts/run_boltz_glyco_sensitivity.sh
venv/bin/python scripts/analyze_boltz_complexes.py --set glyco_sensitivity
venv/bin/python scripts/summarize_glyco_sensitivity.py
venv/bin/python scripts/rank_formal_msa_models.py
bash scripts/run_relax_empty_msa_representatives.sh
bash scripts/run_relax_msa_representatives.sh
venv/bin/python scripts/assess_relaxed_geometry.py
venv/bin/python scripts/analyze_pose_glyco_geometry.py
venv/bin/python scripts/summarize_phase2.py
```

Formal Boltz settings: Boltz-1, seeds 1-15 per route, 3 recycling steps, 200 sampling steps, one diffusion sample, no optional kernels, PDB output, full PAE output. These 30 formal inputs explicitly use empty MSAs.

The sensitivity comparison uses matched 50-step settings and seeds 1-3 for both empty-MSA and ColabFold-MSA conditions. Raw MSA depths and hashes are recorded in `data/msa/msa_manifest.json`; Boltz parsing used the CLI default cap of 8192 sequences per chain. This analysis demonstrates sensitivity to MSA choice but is not a 15-seed replacement formal run.

The completed formal MSA comparison uses seeds 1-15, 3 recycling steps and 200 sampling steps for both Native and Immunogen. The N74 glycosylation sensitivity uses one covalently linked NAG at local antigen residue 48 and the same 200-step settings for seeds 1-3; it is a minimal core-sugar test, not a full human glycoform model.

The Native-Guided route is implemented as Phase1-informed post-ranking of the Native-Blind ensemble because a fair two-chain antibody pocket restraint was not established. It is not described as directly restrained prediction.

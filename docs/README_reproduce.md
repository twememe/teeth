# Reproduction guide

## Minimum system requirements

Use Linux, Python 3.10.x, CUDA 12.4, torch 2.6.0+cu124, an RTX A6000-class GPU, 24 GiB VRAM, 32 GiB RAM, HMMER and at least 20 GiB free disk beyond model assets. Install `requirements-full.txt` or `environment.yml`.

## Quick result reproduction

Run `bash run.sh`. It reconstructs root `results.csv` from the frozen formal-MSA ranking and attached PDBs without GPU work or model downloads.

## Model download and verification

Run `python3 scripts/fetch_model_assets.py --list`, then `python3 scripts/fetch_model_assets.py`. Boltz uses the verified cache directly. ImmuneBuilder 1.2 does not accept this project cache as a model-directory argument; run `python3 scripts/fetch_model_assets.py --install-immunebuilder` to atomically install verified weights into the package's actual `trained_model` directory. Confirm all runtime locations with `python3 scripts/fetch_model_assets.py --verify-only`.

## Full formal GPU reproduction

Run `bash run_full.sh --preflight-only` and resolve every reported item. Then run `bash run_full.sh`. The entry creates a new timestamped directory under `work/`, exports it as `PHASE2_OUTPUT_ROOT`, records every command and exit status, and stops immediately on failure. The formal path calls `run_imgt_numbering.py`, `step11_antibody_structure.py`, `generate_boltz_inputs.py`, `run_boltz_formal_msa.sh`, `analyze_boltz_complexes.py --set formal_msa`, `run_prodigy_formal.py --set formal_msa`, `rank_formal_msa_models.py`, `summarize_formal_msa.py` and `build_recomputed_results.py`.

## Inputs and MSA

IL-24 Native 27–181, Immunogen 27–160, full reference 1–181 and IA6-13-8 VH/VL are frozen under `data/`. Six ColabFold-MSA files and their hashes are recorded by `data/msa/msa_manifest.json`. The MSA is prediction input, not a hidden evaluation set. Formal ColabFold-MSA ensembles use seeds 1–15 per route, three recycling steps, 200 sampling steps and one diffusion sample.

## Supplemental analyses

MSA sensitivity uses matched seeds 1–3 and 50 sampling steps. N74-NAG sensitivity uses seeds 1–3 and 200 sampling steps with one core GlcNAc. These are supplemental uncertainty analyses, not final candidate-generation steps. To rerun them, set `PHASE2_OUTPUT_ROOT` to a new work directory, generate Boltz inputs, then run the corresponding sensitivity shell script, `analyze_boltz_complexes.py` set and summarizer. Restrained relaxation uses `run_relax_msa_representatives.sh` followed by `assess_relaxed_geometry.py`; it is supplemental geometry QC.

## Contact PRODIGY and final results

Contact analysis uses 4.5 Å and 5.0 Å thresholds and project residue numbering. PRODIGY 2.4.0 supplies relative interface-energy ranking only, not experimental KD. Ranking copies ten isolated candidate PDBs and writes the run directory's `results.csv`; no run automatically overwrites the frozen root table.

## Runtime errors and limits

Preflight verifies exact package versions, CUDA, all six model assets, actual ImmuneBuilder runtime weights, MSA hashes, inputs, HMMER, scripts, disk space, path portability and absence of temporary artifacts. Missing items return nonzero and create no output directory. Raw structures may contain clashes; candidate contacts are not experimentally confirmed epitopes. Phase-1 hotspots are priors. External 6DF3 and ELISA are posterior analyses. The recorded 30-job formal-MSA run took approximately 39 minutes on an RTX A6000, excluding preparation and downstream analysis.

# IL-24 antibody structure release package

This package reproduces and audits the Phase 2 computational ranking of IA6-13-8 antibody complexes with IL-24. It provides a deterministic offline result rebuild and a separate GPU recomputation path. The project did not train or fine-tune Boltz-1 or ImmuneBuilder; its contribution is the frozen-input ensemble workflow, interface analysis, sensitivity analyses and auditable ranking.

## Quick result reconstruction

Run `bash run.sh` from the release root. It reads only `results/model_ranking_formal_msa.csv` and the attached candidate PDB files, then writes UTF-8 `results.csv`. It does not download models, use a GPU or modify PDBs and analysis tables. Repeated runs are byte-identical.

The output contains 10 candidates with fixed field order. `track` is `赛道一：AI大分子与多肽药物设计`; `task` is `IL-24 antibody structure prediction`. Numeric missing values are empty. `ipTM` and `complex_pLDDT` are model metrics. `PRODIGY_dG_kcal_mol` is in kcal/mol and is a relative modeling value, not experimental KD. `hotspot_coverage` is a fraction, not proof of an epitope. Ranking is by formal-MSA confidence score and ipTM with recorded confidence and consensus ranks. The complete two-route ranking remains in `results/model_ranking_formal_msa.csv`; route-split result files are intentionally omitted to avoid ambiguous empty outputs.

## Full formal GPU recomputation

Install the exact environment in `requirements-full.txt` or `environment.yml`, install HMMER, obtain models, and run:

```bash
python3 scripts/fetch_model_assets.py --list
python3 scripts/fetch_model_assets.py
python3 scripts/fetch_model_assets.py --install-immunebuilder
python3 scripts/fetch_model_assets.py --verify-only
bash run_full.sh --preflight-only
bash run_full.sh
```

The full entry runs the actual formal workflow: IMGT numbering, ImmuneBuilder ABodyBuilder2, Boltz input generation, the two 15-seed formal ColabFold-MSA ensembles, contact analysis, PRODIGY, ranking, formal summary and isolated `results.csv` generation. It writes only to a new `work/full_recompute_<UTC timestamp>/` directory and stops on the first failure. It never promotes or overwrites frozen root results.

Boltz reads its checkpoint and CCD from `models/boltz/` by default, or from `$PHASE2_MODEL_CACHE/boltz/`. ImmuneBuilder reads its four weights from the installed package's `ImmuneBuilder/trained_model/` directory. `--install-immunebuilder` copies already hash-verified cache files to that actual runtime directory and verifies the installed copies atomically. Set `PHASE2_HMMER_BIN` to the directory containing `hmmscan` when HMMER is not on `PATH`.

Formal settings are seeds 1–15 per Native-Blind and Immunogen-Blind route, three recycling steps, 200 sampling steps and one diffusion sample. MSA sensitivity (matched seeds 1–3, 50 sampling steps) and N74 single-NAG sensitivity (seeds 1–3, 200 sampling steps) are supplemental analyses and do not generate final candidates. Restrained relaxation is supplemental geometry QC. Reproduction commands are in `docs/README_reproduce.md`.

## Inputs outputs and resources

Frozen inputs are under `data/`; every file is recorded in `docs/DATA_PROVENANCE.csv`. Model and software provenance are in `docs/MODEL_PROVENANCE.csv` and `docs/SOFTWARE_VERSIONS.csv`. Structure conventions are in `docs/STRUCTURE_FILES.md`. The compact release omits multi-GB model assets but records exact sizes, hashes and sources in `models/MODEL_ASSETS.csv`.

Minimum full-run environment: Linux, Python 3.10.x, CUDA 12.4, torch 2.6.0+cu124, an NVIDIA RTX A6000-class GPU with at least 24 GiB VRAM, 32 GiB system RAM and at least 20 GiB free disk beyond model assets. The recorded 30 formal-MSA jobs took about 39 minutes on one RTX A6000; preparation and analysis add time. Quick reconstruction requires neither GPU nor model assets.

## Interpretation limits

All PDBs are computational models, not crystal structures, and raw Boltz coordinates may contain severe clashes. Phase-1 hotspots are ranking priors, not epitope truth. The N74 model contains one core GlcNAc, not a complete human glycoform. External 6DF3 and ELISA data are posterior consistency materials and were not training data. 6DF3 represents IL-24–IL-22R1–IL-20R2 and cannot establish an IL-20R1–IL-20R2 complex. Native-Guided denotes Phase-1-informed post-ranking of Native-Blind, not directly restrained prediction.

## Verification

Run `python3 -B validate_release.py`, `bash run.sh`, `python3 -B scripts/fetch_model_assets.py --verify-only`, `bash run_full.sh --preflight-only`, and `sha256sum -c RELEASE_SHA256SUMS.txt`. A missing model or dependency must make full preflight return nonzero before any output directory is created. The optional Notebook is omitted because both documented entry points cover the complete quick and formal workflows.

Source project HEAD: `832af53d02341379637544cac3900f840f763e36`. Source work-tree cleanliness is not asserted because its filesystem was full during review.

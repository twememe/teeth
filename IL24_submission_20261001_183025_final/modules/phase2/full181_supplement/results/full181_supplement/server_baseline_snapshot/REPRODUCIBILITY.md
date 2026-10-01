# Reproducibility record

## Workflow type

No model was trained or fine-tuned for this project. The workflow applies pretrained models, deterministic analysis scripts and ensemble ranking.

## Randomness and parameters

Formal Boltz ensembles use seeds 1-15 per Native and Immunogen route, three recycling steps, 200 sampling steps and one diffusion sample. MSA sensitivity uses matched seeds 1-3 and 50 sampling steps. N74 single-NAG sensitivity uses seeds 1-3 and 200 sampling steps.

## Hardware and time

The recorded environment is Linux, Python 3.10.12, CUDA 12.4, torch 2.6.0+cu124 and NVIDIA RTX A6000. Per-run logs are retained on the server; complete package versions are in `docs/pip_freeze_final.txt`.

## Data separation and leakage

Native and Immunogen constructs are evaluated separately and compared only after prediction. Phase 1B candidate regions are treated as a prior for post-ranking, not prediction ground truth. No hidden competition evaluation set was used. External 6DF3 and ELISA data were introduced for post-hoc mechanistic consistency analysis and were not used to train or fine-tune a model.

## Expected outputs

`bash run.sh` regenerates `results.csv` from the frozen formal-MSA ranking and bundled top-ranked PDB structures. Full GPU recomputation is described separately in `docs/README_reproduce.md` because it is substantially more expensive and requires installing the recorded runtime/model assets.

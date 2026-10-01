# Phase 2 Model Card

## Model and versions

- ImmuneBuilder ABodyBuilder2 1.2 for antibody Fv structure prediction.
- Boltz-1 2.2.1 for IL-24-antibody complex prediction.
- PRODIGY 2.4.0 for relative interface-energy ranking.

## Inputs and outputs

Inputs are the project IL-24 Native/Immunogen FASTA sequences, IA6-13-8 VH/VL sequences, optional ColabFold MSAs, and prediction YAML files. Outputs include PDB structures, interface contacts, confidence metrics, relative PRODIGY scores, ranked representatives, and `results.csv`.

## Reproduction

The completed result list is regenerated with `bash run.sh` (or `python3 predict.py`). Full prediction and analysis commands are documented in `docs/README_reproduce.md`; the compact release includes top-ranked PDB structures and does not include multi-gigabyte inference caches.

## Scope and limitations

This is an inference and ranking workflow, not an experimentally validated epitope or affinity measurement. PRODIGY values are relative model comparisons, not experimental KD. MSA choice, seed-level pose heterogeneity, N74 glycosylation, and incomplete receptor data affect uncertainty. 6DF3 represents IL-24-IL-22R1-IL-20R2 and does not establish the IL-20R1-IL-20R2 complex.

## Innovation contribution

The project contribution is the paired Native/Immunogen ensemble workflow, Phase 1B candidate-region re-evaluation, MSA and N74 sensitivity analysis, 6DF3 receptor-interface mapping, and auditable multi-criterion posterior ranking.

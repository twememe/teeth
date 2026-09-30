# Model card

- ImmuneBuilder ABodyBuilder2 1.2: antibody Fv structure prediction.
- Boltz-1 2.2.1: IL-24–antibody complex prediction.
- PRODIGY 2.4.0: relative interface-energy ranking.

Inputs are IL-24 Native 27–181, Immunogen 27–160, full reference sequence, IA6-13-8 VH/VL, optional ColabFold MSA and prediction configuration. Outputs are PDB structures, contact tables, confidence metrics, quality/QC records and ranked CSV files. Formal Boltz runs use seeds 1–15, 3 recycling steps, 200 sampling steps and one diffusion sample. MSA sensitivity uses matched seeds 1–3 and 50 sampling steps. N74 sensitivity uses one core NAG and seeds 1–3 at 200 sampling steps.

The recorded hardware is NVIDIA RTX A6000, Linux, Python 3.10.12, CUDA 12.4 and torch 2.6.0+cu124. The project did not train or fine-tune the foundation models; it uses open pretrained models, deterministic analysis scripts and ensemble ranking. Native-Guided is Phase-1-informed post-ranking of Native-Blind, not a directly restrained prediction.

Scope is hypothesis generation. ipTM and complex pLDDT are model metrics; PRODIGY values are relative modeling values, not experimental KD; hotspot coverage is not epitope proof. MSA choice, seed-level pose heterogeneity, N74 glycosylation simplification and missing receptor-complex data are uncertainty sources. All PDBs are computational models and may contain clashes. Licenses and sources are in `licenses/`, `docs/MODEL_PROVENANCE.csv` and `THIRD_PARTY_NOTICES.md`.

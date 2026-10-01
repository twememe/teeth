# Third-party software and data

The project uses pretrained/open-source software rather than project-specific model training.

- Boltz-1 2.2.1: biomolecular structure prediction; source and license are retained under `vendor/boltz/` in the working environment.
- ImmuneBuilder 1.2 / ABodyBuilder2: antibody Fv structure prediction.
- ANARCI 2026.2.13.2: antibody numbering.
- PRODIGY 2.4.0: protein-interface affinity estimation used only for relative ranking.
- OpenMM 8.1.1: restrained minimization.
- ColabFold MSA service: MSA generation; raw depths and hashes are in `data/msa/msa_manifest.json`.
- RCSB PDB 6DF3: public experimental structure, downloaded from https://www.rcsb.org/structure/6DF3.

Exact installed dependencies are recorded in `docs/pip_freeze_final.txt`. Users must review and comply with each upstream license before redistribution. The supplied papers retain their publishers' licensing terms and are archived as reference material, not relicensed by this project.

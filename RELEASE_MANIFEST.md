# Release manifest

This compact release contains source code, frozen inputs and MSAs, provenance records, two model licenses, 73 indexed successful-run logs, a 30-model two-route formal-MSA ranking, ten standardized candidates and their PDBs, supplemental analysis tables, entry points and whole-package hashes. It excludes environments, Git metadata, caches, temporary files, complete raw inference trees, unverifiable paper PDFs and multi-GB model weights.

`run.sh` rebuilds the frozen final table only. `run_full.sh` executes the actual formal Boltz workflow and writes solely to a new child of `work/`. Supplemental MSA, N74-NAG and restrained-relaxation scripts also honor `PHASE2_OUTPUT_ROOT`. Model assets are not claimed present unless runtime verification succeeds; exact metadata is in `models/MODEL_ASSETS.csv`.

Source project HEAD is `832af53d02341379637544cac3900f840f763e36`. Source work-tree cleanliness was not asserted because its filesystem was full during review.

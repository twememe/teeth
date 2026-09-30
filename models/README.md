# Model assets

The compact release does not bundle large weights. `MODEL_ASSETS.csv` records six verified assets with byte sizes, SHA-256 and sources. `scripts/fetch_model_assets.py` streams into `.partial` files, validates size and hash, then atomically renames. Boltz uses `models/boltz/` or `$PHASE2_MODEL_CACHE/boltz/` directly.

ImmuneBuilder ABodyBuilder2 1.2 reads weights from the installed package's `ImmuneBuilder/trained_model/` directory, not from the release cache. Download to the cache first, then use `--install-immunebuilder`; this copies only verified files to the runtime directory and verifies them again. `--verify-only` checks the paths the programs actually read. Quick result reconstruction needs no model assets; full GPU recomputation requires all six.

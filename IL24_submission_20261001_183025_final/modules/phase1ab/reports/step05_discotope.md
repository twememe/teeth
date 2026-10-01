# Step 5 — DiscoTope-3.0

## Status

`MISSING_EXTERNAL_DEPENDENCY_FAILURE`

`discotope_available = false`

## What was attempted

Two formal installation attempts were performed against the official DiscoTope-3.0 repository at commit `35d9f2e55f97eaba2a7acefbc394db58fb9670bc`, with the preregistered AlphaFold structure mode reserved for inference.

- Attempt 1 reached the command time ceiling during official CUDA dependency retrieval. The partial dependency transaction was not treated as a usable environment.
- Attempt 2 reused the same dependency plan and cache, but WSL network/DNS failure prevented retrieval from `pypi.org` / `files.pythonhosted.org`. `torch_geometric`, build-time `setuptools`, and torch were not installed; official inference could not start.

## Preserved evidence

- `logs/install_dt3_attempt1.log`
- `logs/install_dt3_attempt2.log`

## Preregistered handling

The two-attempt ceiling is exhausted. No third installation was attempted, no alternative model was substituted, and no pseudo 3D epitope score was constructed. Step 5 is missing but is not a Phase 1 stop condition. Step 7 therefore uses the explicitly authorized available-evidence consensus across BepiPred, RSA, and conservation only.

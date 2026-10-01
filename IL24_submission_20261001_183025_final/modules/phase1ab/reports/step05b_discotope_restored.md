# Step 05B — Restored DiscoTope-3.0 evidence

## Status

- Status: **COMPLETE** (`discotope_available = true`).
- Restoration authority: post-freeze Phase 1B amendment; Phase 1A remains immutable.
- Historical experimental regions were not accessed.

## Official method and installation provenance

- Official repository: https://github.com/DTU/DiscoTope-3.0
- Frozen commit: `35d9f2e55f97eaba2a7acefbc394db58fb9670bc`
- Official README runtime: Python 3.14; actual isolated environment: `il24-dt3-phase1b`, Python 3.14.7.
- B1 preserved the official dependency specifications. Network Range retries completed the large wheels, after which `biotraj` produced an explicit source-build failure because `gcc` was absent.
- B2 added only the compatible Conda GCC/G++ 16.1.0 toolchain. No predictor, model, Python version, or official requirement version was substituted.
- Official loaders successfully loaded 100 XGBoost models, two GAM calibration models, and ESM-IF1 (141,662,150 parameters).

DiscoTope-3.0 uses an ESM inverse-folding protein model to encode the three-dimensional backbone and sequence context, then applies an ensemble of XGBoost residue classifiers. The official length/surface GAM models calibrate scores across proteins. The supplied AlphaFold structure was therefore run with `--struc_type alphafold`.

## Inference

- Device actually used: NVIDIA GeForce RTX 4050 Laptop GPU (CUDA 12.8); no CUDA OOM and no CPU fallback.
- Wall-clock time: 55.83 s.
- Observed GPU memory: 1587–3954 MiB (peak delta 2367 MiB across 56 samples).
- Input gate: one model, Chain A only, residues 1–181, complete backbone, exact reference sequence match.
- Official raw output: exactly 181 rows; no non-finite continuous values.

## Residue scores

- Official uncalibrated `DiscoTope-3.0_score` is preserved as `DiscoTope_raw` and is the continuous value preregistered for Phase 1B percentile ranking.
- Official `calibrated_score` is separately preserved as `DiscoTope_calibrated`; the official epitope classification uses its default threshold 0.90.
- `DiscoTope_raw` range: 0.00235–0.57638.
- Calibrated score range: -1.39816–2.49046.
- Official positive residues: 38/181; contiguous positive segments: 30-30, 33-33, 46-46, 49-49, 54-54, 58-58, 60-60, 63-64, 68-68, 73-74, 76-76, 93-93, 97-98, 102-102, 105-106, 110-110, 113-113, 117-117, 120-120, 124-124, 127-127, 130-135, 137-137, 139-139, 142-143, 150-150, 154-154, 173-173.
- Ten highest raw-score residues: R73 (0.57638), K130 (0.53788), D131 (0.50583), Q143 (0.49313), S33 (0.48055), R110 (0.47963), R150 (0.47883), N117 (0.46785), W49 (0.46289), N93 (0.44856).
- Highest nonredundant 15-aa raw-score windows: 129-143 (mean 0.36944), 121-135 (mean 0.30823), 97-111 (mean 0.28260).

These score summaries are derived only from the restored predictor output and do not use historical assay intervals. Formal multi-evidence candidate selection is deferred to Steps 7B–8B.

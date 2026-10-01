# IL-24 Phase 1B Pre-unblinding Report

## 1. Why Phase 1B was performed

The preregistered DiscoTope-3.0 structural-epitope branch was missing in Phase 1A because the allowed installation attempts ended in external dependency failure. Before any historical experimental unblinding, the user explicitly authorized a post-freeze restoration of that missing branch. This is a protocol amendment, not a retroactive claim that the restoration was part of the original successful run.

## 2. Frozen Phase 1A

Phase 1A outputs remained read-only. The original freeze manifest was replayed after Phase 1B: 68/68 entries matched, zero mismatches, manifest SHA-256 `23c0330bf6bba11f55727689cfe27fbcdb3a4b125829314aa046614166490a5a`. Historical experimental regions were not accessed.

## 3. DiscoTope-3.0

Official DiscoTope-3.0 commit `35d9f2e55f97eaba2a7acefbc394db58fb9670bc` ran in `il24-dt3-phase1b` (Python 3.14.7, Torch 2.10.0/CUDA 12.8). B2 added only GCC/G++ to compile the official `biotraj` dependency; predictor code, weights, and requirement versions were unchanged. Official ESM-IF1, 100 XGBoost models, and two GAM calibrators loaded successfully. AlphaFold mode used the RTX 4050 GPU without OOM or CPU fallback; wall time was 55.83 s.

DiscoTope is a structure-aware B-cell conformational epitope predictor. Unlike a sequence-only model, it encodes the folded backbone so residues far apart in the linear chain can jointly contribute to an antibody-accessible 3D surface patch. RSA asks whether a residue is exposed; DiscoTope asks whether that residue and its 3D neighborhood resemble an epitope.

## 4. 3D epitope results

The official output mapped exactly to positions 1-181 with no non-finite values. `DiscoTope_raw` is the official uncalibrated continuous score used for percentile ranking; official calibrated scores and the default 0.90 classification were separately preserved. Raw range: 0.00235-0.57638; official positive calls: 38/181.

## 5. BepiPred vs DiscoTope

Raw Pearson r = 0.6271; raw Spearman rho = 0.6426. The moderate agreement indicates both convergence and complementary information rather than redundancy.

## 6. Full evidence integration

All features use the Phase 1A average-tie percentile transform on [0,1]. `E=(B+D)/2`, then the primary score is `FinalScore_Phase1B=(E+R+C)/3`, giving one-third weight each to epitope consensus, surface accessibility, and evolution. `(B+D+R+C)/4` is sensitivity-only. pLDDT is a confidence flag and DSSP is annotation; neither changes the score.

## 7. Phase 1B candidate regions

All 15-, 20-, and 25-aa windows were freshly recomputed from FinalScore Phase 1B. The unchanged family rule is overlap coefficient `intersection / shorter length > 0.5`; no single scale or unique Top1 was selected.

| Candidate | Region | FinalB | BepiPred | DiscoTope | Epitope consensus | RSA | Conservation | median pLDDT | pLDDT>=70 | Scale | LOO | Bootstrap median [95% CI] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| C1-B | 52-66 | 0.6073 | 0.2050 | 0.2634 | 0.6957 | 0.4588 | 0.8383 | 92.19 | 1.00 | 3/3 (15;20;25) | 4/4 | 0.8520 [0.7977, 0.8982] |
| C2-B | 126-140 | 0.5855 | 0.1684 | 0.3504 | 0.7098 | 0.4913 | 0.7845 | 90.25 | 1.00 | 3/3 (15;20;25) | 3/4 | 0.7880 [0.7238, 0.8652] |
| C3-B | 60-74 | 0.5693 | 0.1478 | 0.2696 | 0.5937 | 0.4885 | 0.8169 | 96.06 | 1.00 | 3/3 (15;20;25) | 3/4 | 0.8232 [0.7672, 0.8775] |
| C4-B | 143-157 | 0.5556 | 0.1582 | 0.2667 | 0.6122 | 0.3726 | 0.8171 | 98.19 | 1.00 | 3/3 (15;20;25) | 2/4 | 0.8775 [0.8232, 0.9439] |
| C5-B | 117-131 | 0.5468 | 0.1798 | 0.2689 | 0.6472 | 0.3715 | 0.8078 | 93.81 | 1.00 | 3/3 (15;20;25) | 2/4 | 0.8477 [0.7496, 0.9413] |

## 8. Robustness

Grouped leave-one-evidence-out results are reported as matched family rank, matched mean, and overlap. Scale stability uses region-level matching across all three fixed scales. Conservation uncertainty reuses the frozen Step 6 alignment with exactly 100 bootstrap repetitions.

| Candidate | Full | -Epitope | -Surface | -Conservation |
|---|---:|---:|---:|---:|
| C1-B | rank 1; mean 0.6073; ov 1.00 | rank 3; mean 0.5674; ov 0.73 | rank 1; mean 0.6258; ov 1.00 | rank 2; mean 0.6413; ov 0.80 |
| C2-B | rank 2; mean 0.5855; ov 1.00 | rank 8; mean 0.5261; ov 0.93 | rank 2; mean 0.5836; ov 0.73 | rank 1; mean 0.6617; ov 0.80 |
| C3-B | rank 3; mean 0.5693; ov 1.00 | rank 3; mean 0.5674; ov 0.73 | rank 5; mean 0.5578; ov 1.00 | rank 6; mean 0.6006; ov 0.80 |
| C4-B | rank 4; mean 0.5556; ov 1.00 | rank 7; mean 0.5272; ov 1.00 | rank 3; mean 0.5827; ov 0.93 | rank 8; mean 0.5609; ov 1.00 |
| C5-B | rank 5; mean 0.5468; ov 1.00 | rank 10; mean 0.5085; ov 1.00 | rank 4; mean 0.5591; ov 0.80 | rank 9; mean 0.5120; ov 0.87 |

Epitope-model sensitivity:

- BepiPred-only Top 5: 52-66 (rank 1, primary match C1-B, ov 1.00), 123-137 (rank 2, primary match C2-B, ov 0.80), 60-74 (rank 3, primary match C3-B, ov 1.00), 145-159 (rank 4, primary match C4-B, ov 0.87), 24-38 (rank 5, primary match NO_TOP5_MATCH, ov 0.00)
- DiscoTope-only Top 5: 129-143 (rank 1, primary match C2-B, ov 0.80), 52-66 (rank 2, primary match C1-B, ov 1.00), 60-74 (rank 3, primary match C3-B, ov 1.00), 142-156 (rank 4, primary match C4-B, ov 0.93), 121-135 (rank 5, primary match C5-B, ov 0.73)

The BepiPred-only residue score reproduces Phase 1A to numerical precision. The DiscoTope-only and combined rankings were computed independently; no weights were tuned to preserve a Phase 1A hotspot.

## 9. Phase 1A vs Phase 1B

| Phase 1A | Phase 1B match | Overlap | Rank A -> B | 3D delta (A region) |
|---|---|---:|---:|---:|
| C1-A 52-66 | C1-B 52-66 | 1.00 | 1 -> 1 | -0.0233 |
| C2-A 123-137 | C2-B 126-140 | 0.80 | 2 -> 2 | +0.0033 |
| C3-A 60-74 | C3-B 60-74 | 1.00 | 3 -> 3 | +0.0054 |
| C4-A 145-159 | C4-B 143-157 | 0.87 | 4 -> 4 | -0.0117 |
| C5-A 24-38 | NO_TOP5_MATCH | 0.00 | 5 -> not Top 5 | -0.0245 |

## 10. Interpretation

Regions with high mean BepiPred and DiscoTope percentile have convergent sequence-based and structure-based epitope support. Regions reduced by Phase 1B had strong sequence evidence but weaker conformational support; newly promoted regions gained structure-supported priority. These are all acceptable pre-unblinding outcomes.

- C1-B (52-66) and C3-B (60-74) may be interpreted as a broader biological hotspot super-family; formal families remain separate.
- C2-B (126-140) and C4-B (143-157) may be interpreted as a broader biological hotspot super-family; formal families remain separate.
- C2-B (126-140) and C5-B (117-131) may be interpreted as a broader biological hotspot super-family; formal families remain separate.

Super-family language is interpretive only: it does not merge candidates or alter the fixed clustering threshold.

## 11. Limitations

- These are not experimentally validated epitopes.
- These are not receptor interfaces.
- These are not confirmed neutralizing epitopes.
- These are not periodontitis-specific functional interfaces.
- AlphaFold confidence and computational predictor scores do not establish biological function.

## 12. Ready for unblinding

Phase 1B pre-unblinding analysis frozen. Historical experimental intervals remain blinded; Step 9 was not performed.

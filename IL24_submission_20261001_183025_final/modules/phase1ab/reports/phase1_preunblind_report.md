# IL-24 Phase 1 Pre-unblinding Report

## 1. Research question

Prioritize continuous regions on the frozen 181-aa IL-24 project sequence using only preregistered evidence available before historical experimental unblinding.

## 2. Frozen input

The coordinate system is project positions 1–181 from `inputs/reference_181.fasta`. The AlphaFold structure and UniProt mapping input retain their Step 0 manifest hashes. No historical experimental interval or excluded image was opened or used.

## 3. Sequence QC

The frozen reference, PDB Chain A, and all downstream structural residue tables map continuously across 181 project residues. The formal Step 1 output was integrity-audited but not rerun.

## 4. Structure QC

The AlphaFold Chain A structure contains 181 mapped residues. pLDDT is used only as confidence metadata and never contributes to FinalScore.

## 5. BepiPred-3.0

Official standalone BepiPred-3.0 `vt_pred` uses ESM-2 protein-language-model representations learned from amino-acid sequence context, analogous to a language model learning contextual word patterns. The formal run used RTX 4050 CUDA without CPU fallback and completed in 84.504 s. The output contains 181 residue probabilities; range 0.018547–0.313823. Continuous probabilities, not only threshold calls, enter the consensus.

## 6. Surface accessibility

FreeSASA RSA provides independent physical surface-accessibility evidence. DSSP secondary structure is annotation only and does not add or subtract score.

## 7. DiscoTope-3.0

`discotope_available = false` and status is `MISSING_EXTERNAL_DEPENDENCY_FAILURE`. Both preregistered official installation attempts were preserved; Attempt 2 failed during WSL DNS/PyPI dependency retrieval. No third attempt, substitute model, or artificial 3D epitope score was used. This preregistered structural-epitope evidence branch is therefore missing.

## 8. Evolutionary conservation

MAFFT-aligned mammalian orthologs were scored with normalized Shannon conservation. Candidate-level uncertainty uses exactly 100 ortholog bootstrap resamples.

## 9. Consensus methodology

The primary analysis is a transparent available-evidence consensus. BepiPred probability, RSA, and conservation were independently average-tie percentile ranked to [0,1], then combined as `FinalScore = (B + R + C) / 3`. pLDDT is a structure-confidence flag; DSSP is annotation.

## 10. Candidate regions

All 15-, 20-, and 25-aa windows were retained. Score-ordered interval families use overlap coefficient `intersection / shorter length > 0.5`; no single window size or unique Top1 was selected.

| Candidate | Region | Length | Final | BepiPred | RSA | Conservation | median pLDDT | pLDDT >=70 | Structural confidence | Scale | LOO | Bootstrap median [95% CI] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| C1 | 52–66 | 15 | 0.6306 | 0.2050 | 0.4588 | 0.8383 | 92.19 | 1.00 | HIGH | 3/3 (15;20;25) | 4/4 | 0.8520 [0.7977, 0.8982] |
| C2 | 123–137 | 15 | 0.5819 | 0.1860 | 0.4758 | 0.7723 | 90.25 | 1.00 | HIGH | 3/3 (15;20;25) | 2/4 | 0.7926 [0.7374, 0.8652] |
| C3 | 60–74 | 15 | 0.5638 | 0.1478 | 0.4885 | 0.8169 | 96.06 | 1.00 | HIGH | 3/3 (15;20;25) | 3/4 | 0.8232 [0.7672, 0.8775] |
| C4 | 145–159 | 15 | 0.5624 | 0.1657 | 0.3770 | 0.8143 | 97.62 | 1.00 | HIGH | 1/3 (15) | 2/4 | 0.8775 [0.8232, 0.9439] |
| C5 | 24–38 | 15 | 0.5619 | 0.1513 | 0.4049 | 0.8343 | 79.62 | 0.80 | HIGH | 3/3 (15;20;25) | 3/4 | 0.9069 [0.8520, 0.9439] |

## 11. Robustness analysis

Each Full candidate was region-matched against independently recomputed Full, -Epitope, -Surface, and -Conservation families. The stability matrix contains 20 candidate-setting records with rank, mean score, and region overlap. Scale stability is defined by overlap with top-five scale-specific families rather than exact endpoints.

Cells are `rank / mean / overlap` for the matched region family.

| Candidate | Full | -Epitope | -Surface | -Conservation |
|---|---:|---:|---:|---:|
| C1 | 1 / 0.6306 / 1.00 | 3 / 0.5674 / 0.73 | 1 / 0.6607 / 1.00 | 1 / 0.6774 / 0.80 |
| C2 | 2 / 0.5819 / 1.00 | 8 / 0.5261 / 0.73 | 6 / 0.5562 / 0.87 | 2 / 0.6567 / 1.00 |
| C3 | 3 / 0.5638 / 1.00 | 3 / 0.5674 / 0.73 | 8 / 0.5496 / 1.00 | 4 / 0.6200 / 0.80 |
| C4 | 4 / 0.5624 / 1.00 | 7 / 0.5272 / 0.87 | 4 / 0.5871 / 1.00 | 6 / 0.5743 / 1.00 |
| C5 | 5 / 0.5619 / 1.00 | 4 / 0.5482 / 1.00 | 3 / 0.5983 / 0.93 | 10 / 0.5102 / 0.73 |

## 12. Limitations

The DiscoTope structural-epitope branch is missing, BepiPred thresholds are not protein-specific calibration, AlphaFold confidence is not functional evidence, and conservation depends on the curated ortholog set. These results are multi-evidence computational candidate regions—not experimentally validated epitopes, receptor interfaces, or confirmed neutralizing epitopes.

Phase 1 remains pre-unblinded. Step 9 was not performed.

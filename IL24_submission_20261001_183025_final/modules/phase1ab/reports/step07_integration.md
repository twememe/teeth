# Step 7 — Multi-Evidence Integration

## Status

`discotope_available = false`. The preregistered DiscoTope-3.0 structural-epitope branch is missing because both allowed official installation attempts ended in external dependency retrieval failure. No substitute model or fabricated 3D epitope score was introduced.

## Available-evidence consensus

Each available continuous feature was independently transformed with average-tie percentile ranks scaled to [0, 1]: BepiPred-3.0 probability (B), FreeSASA RSA (R), and normalized Shannon conservation (C).

`FinalScore = (B + R + C) / 3`

pLDDT is retained only as a structure-confidence flag. DSSP is retained only as annotation; neither contributes to FinalScore.

## Validation and summary

- Rows: 181; project positions: 1–181.
- FinalScore range: 0.129630–0.813889; mean 0.500000.
- BepiPred–RSA Pearson r: 0.4812.
- BepiPred–Conservation Pearson r: -0.3243.
- RSA–Conservation Pearson r: -0.5178.

These are multi-evidence computational prioritization scores, not experimentally validated epitopes, receptor interfaces, or confirmed neutralizing regions.

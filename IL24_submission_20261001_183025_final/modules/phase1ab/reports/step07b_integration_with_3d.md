# Step 7B — Complete integration with restored 3D evidence

## Primary formula

Each raw feature uses the frozen Phase 1A average-tie percentile transform on [0,1]. `DiscoTope_percentile` is computed from the official uncalibrated `DiscoTope-3.0_score` preserved as `DiscoTope_raw`.

`EpitopeConsensus = (BepiPred_percentile + DiscoTope_percentile) / 2`

`FinalScore_Phase1B = (EpitopeConsensus + RSA_percentile + Conservation_percentile) / 3`

The direct four-feature mean is retained only as `FinalScore_4way_sensitivity` and is not used for primary candidate ranking. pLDDT is a confidence flag only; DSSP is annotation only.

## Validation

- Rows: 181; reference/numbering exact.
- BepiPred-only reproduction of frozen Phase 1A: maximum absolute error 8.327e-17.
- Primary Phase 1B score range: 0.129630–0.822222.

## BepiPred versus DiscoTope

- Raw Pearson r: 0.6271.
- Raw Spearman rho: 0.6426.
- Percentile Pearson r: 0.6426.
- Percentile Spearman rho: 0.6426.
- Extreme cross-model conflict residues (one >=0.8 while the other <=0.2): 1/181.

Low or moderate correlation is not a failure: the sequence model and conformational structure model interrogate complementary epitope properties.

## Phase 1A versus Phase 1B residue changes

- Delta range: -0.114815 to 0.100926; mean 0.000000.
- Largest positive deltas: L8 (+0.1009), S113 (+0.0907), N74 (+0.0824), S33 (+0.0741), P9 (+0.0704), I137 (+0.0704), L72 (+0.0694), F90 (+0.0694), F145 (+0.0694), L135 (+0.0676).
- Largest negative deltas: I103 (-0.1148), E28 (-0.1120), Y99 (-0.1065), S31 (-0.1028), Q125 (-0.1019), L126 (-0.0944), Q56 (-0.0898), P128 (-0.0731), Q27 (-0.0676), G77 (-0.0667).

No historical experimental interval was accessed. These are computational score comparisons only.

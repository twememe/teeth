# Baseline audit

Historical documentation records 15 seeds per construct and condition (60 models). The current selected inputs contain 60 saved per-model summaries (empty=30, formal_msa=30), 4459 saved contact rows, and 60 PDB files representing 60 unique models. The optional server CPU audit verifies 60 accessible original PDB/confidence/PAE bundles. Missing per-model rows remain NA; historical reported completion is separate from current availability.

Original per-condition table directories take precedence over flattened release fallbacks. When a flattened summary is used, its formal-MSA identity is checked against saved ranking metrics. results.csv is only a selected candidate export and never substitutes for an ensemble; an empty results_immunogen.csv does not imply Immunogen was not run.

Original ranking: confidence descending; consensus mean pairwise epitope Jaccard descending then confidence; informed frozen-prior coverage descending then confidence. The older empty-MSA ranking did not contain consensus ranking. Full181 applies the formal-MSA consensus rule separately to each condition, with numeric seed ascending as the fixed final tie break. Prior coverage denominator is the entire frozen prior, not the construct intersection.

Consensus uses >=50% support: 8/15. Legacy Jaccard assigns 1 to two empty epitopes; the exports separately mark both_empty and exclude it from nonempty-union summaries. Incomplete Full181 consensus is provisional only; absent models never count as no-contact models.

Frozen prior CSV status: available. The historical generator lists ranges, but this is not substituted for the missing frozen file. Saved old hotspot scores remain labelled historical. No hash was computed.

The empty-MSA shared 27–160 consensus Jaccard is 1/31=0.032258 from the frozen sets, consistent with formal_msa_comparison.csv; Phase2_summary.md gives 0.031, which is not the shared-region recomputation. Geometry uses only actual supplied CPU audit observations; no subset is extrapolated to 15 models. No old GPU or inference task was launched.

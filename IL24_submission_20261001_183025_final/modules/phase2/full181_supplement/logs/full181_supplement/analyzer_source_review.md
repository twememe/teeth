# Bounded source review — 2026-09-21

One blocking issue remains **if the two Full181 branches run concurrently**:
`run_seed1_downstream()` calls `refresh_manifest()`, which directly truncates and
rewrites the shared `results/full181_supplement/task_manifest.csv`
(`run_full181_supplement.py`, lines 212–235). Each condition-filtered analyzer
still reads that shared file (`analyze_full181_supplement.py`, line 257). The
other branch can rewrite it during this read, producing a partial/corrupt or
empty manifest and falsely blocking a valid formal seed 1. Separate
`--output` directories protect analysis products but do not protect this shared
input. Before enabling concurrent branches, publish the manifest atomically
and coordinate concurrent refreshes, or serialize the branches. No source was
modified in this review.

No other blocking issue was found in the requested scope. Explicit PDB-to-local
and project mappings are used throughout; missing/failed observations export NA
while completed zero-contact models export empty sets and zero counts. Contact
distances use unrounded residue-pair minima and strict `<4.5` / `<5.0` cutoffs;
only antigen A versus H/L protein heavy atoms are considered. CDR/IMGT lookup
uses verified antibody sequence positions. The runner's condition/seed filters
match the analyzer CLI, and per-branch seed-1 outputs are isolated under their
respective task directories.

Evidence read, without rerunning it: the existing
`interface_audit/legacy_raw_subset_interface_qc.csv` contains 13 unique legacy
PDB observations with all four saved-interface agreement flags true. The
previous 9 numbering tests passed. This review ran no tests or models and made
no hash calculations. Full remains 0/30 because the original environment is
unreadable; successful Full end-to-end execution is not claimed.
# Root resolution after review

The shared manifest race was corrected: `refresh_manifest()` now serializes
snapshots with an OS file lock and publishes a uniquely named temporary CSV by
atomic replacement. The shared `write_json()` helper also uses atomic replacement
for task states, so concurrent readers cannot read partly written JSON. The
updated manifest exporter was actually run to retain the current 30 blocked tasks;
no inference, model smoke test, or parallel GPU execution was performed.

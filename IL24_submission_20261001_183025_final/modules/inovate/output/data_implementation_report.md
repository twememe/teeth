# Data implementation report

{
  "raw_counts": {
    "abbench": 215699,
    "abagym": 36541,
    "il6": 1636
  },
  "total_raw": 253876,
  "split_counts": {
    "train": 199999,
    "validation": 20000,
    "test": 33508
  },
  "unique_model_sequences": 178850,
  "rejections": {
    "abagym:sequence_length_above_1022": 338,
    "abagym:threshold_tie_0.5": 30,
    "abbench:training_median_tie": 1
  },
  "abbench_training_median": 7.828436341973866,
  "official_il6_split_changes": 0
}

- Pooled labels are heterogeneous; high numeric DMS score can represent immune escape, not antibody affinity improvement.
- AbBiBench mutation annotations are absent from source; mutation stays empty. Source sequences are already measured variants.
- AbBiBench H and L are concatenated without linker; this introduces an artificial chain boundary.
- AbAgym model input is the longest verified representative mutated PDB chain with observed ATOM residues only, not a complete antibody-antigen complex. Unresolved residues are absent. All equivalent chains mutation sites must verify.
- Exact duplicate sequences are grouped across all sources; sequence identity/antigen/family-disjoint generalization is not assessed. Conflicting assays/labels are retained for feasibility and reported.
- Original IL6 split assignments are retained, with test-over-validation-over-train priority if an identical sequence has conflicting official assignments.
- No oversampling or fabricated records. Split targets apply to eligible rows before AbBiBench median ties; deficits are reported.

Source revisions and hashes: data/download_manifest.json. Rejected rows: data/rejected_records.csv. Smoke files contain 1,000 training and 200 validation rows with all available sources represented.

Independent CSV reload validation passed (AA alphabet, length, labels, hashes, all split isolation, all retained AbAgym mutations, original IL6 split preservation). Seven unit tests passed.

Conflicting sequence-label groups by source: {"abbench": 61030, "abagym": 1203}. Within the same assay: {}. Antigen/assay context is not included in model_sequence, so identical inputs can have contradictory labels; this is a material limitation of pooled sequence-only feasibility metrics.

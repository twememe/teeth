# Full181 interface analysis audit

## Scope and minimal changes

The implementation is isolated in `src/analyze_full181_supplement.py`, substantive regression tests, and this report. Historical scripts, Native155/Immunogen134 inputs, model structures, and summary files remain read-only. Hash computations, smoke tests, plots, GPU predictions, and MSA requests are excluded by the user's current instructions.

The source protocol is `D:/Downloads/IL24_Full181_Codex_supplement_prompt.md`. The intended downstream design reads an explicit model task manifest, keeps all planned records, validates PDB sequences and explicit coordinate maps, calculates antigen A versus Fv H,L contacts and the existing geometry rules, and emits separate confidence/PAE provenance. The ensemble/comparison agent consumes the model and contact tables.

## Initial read-only audit

- Project root: `D:/IL-24(teeth）/phase2_release`. No AGENTS.md was found in the project or its parent tree search.
- Read `MODEL_CARD.md`, `REPRODUCIBILITY.md`, `docs/Phase2_summary.md`, contact, geometry, ranking, consensus, and PRODIGY scripts.
- `analyze_boltz_complexes.py` hardcodes antigen local + 26; its contact table omits IMGT positions and source paths. It selects A versus H,L correctly, but its 5.0 Å query is inclusive without a final strict threshold filter.
- Preserve the geometry rules: adjacent CA >4.5 Å; peptide C-N outside 1.1–1.6 Å; interchain heavy-atom proximity at the legacy 1.8 Å search radius; nonlocal intrachain proximity at 1.8 Å and residue-index separation >1. Any nonzero count marks a warning, never exclusion.
- Preserve confidence ranking descending confidence; consensus ranking descending mean Jaccard then confidence; Phase1-informed ranking descending hotspot coverage then confidence. Resolve ties by ascending numeric seed. Two empty epitopes have legacy Jaccard 1, separately marked as empty and not evidence of binding.
- Frozen prior `results/12_complex_prediction/phase1_prior_deduplicated.csv` is absent in the local release. Phase1 overlap/ranking must stay blocked unless the actual frozen file is supplied. Report intervals are not a substitute.
- Available local old material: 30-row formal-MSA model summary/contact/PRODIGY tables; formal and formal-MSA historical aggregate summaries; five formal-MSA representative PDBs and ten Native empty-MSA PDBs. No local raw confidence JSON/PAE ensemble files were found. These observations do not imply the complete old raw ensembles are available.

## Proposed input manifest contract

Required columns: `construct,condition,seed,source_path,status`. Construct is `Full181`; condition is `empty` or `formal_msa`; seeds are integers 1–15. `source_path` identifies the anticipated PDB even for unfinished tasks and may be absolute or resolved relative to the manifest directory. Optional columns: `confidence_path,pae_path,chain_map_path,chain_index_map,chain_index_map_source,validation_status,failure_reason,pdb_plddt_scale,confidence_plddt_scale`.

Only rows with `status=validated` or `validation_status=passed` are analysis candidates. They are rechecked for required PDB/confidence/PAE files and A/H/L sequences. Missing or invalid rows retain the complete model key with metrics as NA. The default expected plan is exactly 30 unique condition/seed tasks; missing manifest rows are explicitly added as missing tasks and duplicates are rejected.

`chain_index_map` is JSON (for example `{"A":0,"H":1,"L":2}`) only when verified from actual predictor chain metadata; `chain_index_map_source` identifies that evidence. Numeric JSON confidence keys are never mapped from PDB chain ordering alone. Without verified mapping, A-H/A-L metrics remain NA.

## Output contract for the ensemble/comparison agent

Outputs under `results/full181_supplement/`:

- `full181_model_summary.csv`: complete keys; `status`, `completed` (1 only for successfully validated and analyzed), failure reason; PDB/confidence/PAE paths; confidence_score, ptm, iptm, complex_plddt; epitope_project_residues (semicolon list; blank with completed=1 means observed no contact); epitope_residues_4p5, paratope_residues_4p5, cdr_paratope_residues_4p5; hotspot_coverage, hotspot_hits_4p5, hotspot_fraction_in_interface (NA if frozen prior absent); 5.0 Å counterparts; pLDDT scale/provenance; chain pair confidences; geometry_pass and warnings.
- `full181_interface_contacts.csv`: keys; antigen PDB residue number/insertion, local sequence index, project position/resname; antibody chain, PDB residue number/insertion, local sequence index, IMGT position/insertion, CDR/FR; exact min_distance_A and strict contact_4p5/contact_5p0 indicators.
- `full181_paratope.csv`: keys plus one antibody contact residue with CDR/FR/IMGT, min distance, partner antigen project positions.
- `full181_geometry_qc.csv`: keys and the unchanged geometry thresholds/counts, no model removal.
- `terminal_contact_summary.csv`: keys plus 1–26, 27–160, 161–181 contact counts/fractions; common-region epitope sets.

Validation evidence and limitations will be appended after implementation and substantive checks.

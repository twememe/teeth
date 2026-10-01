# Full181 runner source review — 2026-09-21

Reviewed and updated `src/run_full181_supplement.py` only, apart from this review
record. The runner has not been invoked for a prediction, model test, smoke test,
or PRODIGY calculation. No digest was computed. A Python syntax-tree parse was
performed during source review before the final request to avoid further trial
runs; it did not import or execute the runner. Final changes were reviewed as
source. End-to-end execution has not been verified and is not claimed.

## Corrections

- Removed the nonexistent `build_antibody_residue_map` import. Acceptance now
  calls `build_structure_mapping(pdb, "Full181", ROOT, require_complete=True)`.
  It verifies the exact frozen A181/H119/L106 sequences and all three protein
  chains, including explicit residue numbering and insertion-code handling.
- PDB acceptance checks finite ATOM coordinates, one coordinate model, and
  finite PDB pLDDT values on the frozen writer's 0–100 scale. Four required
  confidence JSON fields must be finite numeric values on the 0–1 scale.
- PAE acceptance requires a numeric `pae` array with exact shape 406×406,
  finite entries, and no negative values. Shape, type, range, units, and the
  source-file metadata are saved in each accepted task state.
- Input YAML must match the fixed A/H/L sequences and the exact empty/formal
  MSA condition. Formal-MSA deployment validation replays the retained A/H/L
  paired/unpaired A3M files, checks their current batch query and equal paired
  row counts, applies the original 8192/16384 combination policy, and compares
  every resulting CSV key/sequence row directly with the deployed CSV. Queries,
  aligned column lengths, raw depths, retained paired keys, and combined depths
  must agree. Original acquisition path/mtime fingerprints are not compared
  with server paths. Fresh server metadata is recorded separately in
  `msa/current_input_validation.json` and the per-task provenance. Same-server
  resumption continues to compare that runtime metadata exactly.
- Per-task state validates construct, condition, seed, exact input text,
  sequence identity, checkpoint/CCD/executable/MSA metadata, the complete
  formal command, and accepted artifact metadata. Only Full181 conditions
  empty/formal_msa and seeds 1–15 can be addressed by the runner.
- CLI paths are resolved before recording commands, and Boltz/PRODIGY must
  belong to the currently audited original environment. The fixed Boltz 2.2.1,
  torch 2.6.0+cu124, original checkpoint size, and original CCD size checks remain.
- Seed 1 is a real 200-step formal task. After its artifacts are accepted, the
  runner invokes contact analysis and that seed's original unrelaxed
  A-versus-H,L PRODIGY calculation before advancing to seed 2. A downstream
  execution or parsing failure preserves the raw model and stops progression;
  resumption may redo deterministic downstream work without another prediction.
- PRODIGY records stdout, stderr, command, exit code, timestamps, 25 °C,
  interface contacts, and finite DeltaG/Kd when available. The exact frozen
  exception `No contacts found for selection` is recorded as `no_interface`,
  with blank energy/Kd and zero observed contacts; it does not cause selective
  exclusion or an extra prediction. Other parsing/execution failures stay failed.
- Every `Path.read_text()` in this runner now explicitly uses UTF-8. Text log
  writes explicitly use UTF-8. Monitoring failures no longer orphan a running
  formal prediction merely because one GPU status sample failed.

## Analyzer integration contract

At review time `src/analyze_full181_supplement.py` had not yet been written by
the independent analysis agent. The runner currently invokes it with no CLI
arguments from the release root. It must analyze the accepted raw Full outputs
and write these project-level tables under `results/full181_supplement/`:

- `full181_model_summary.csv`: exactly one Full181 row for the requested
  condition/seed 1, including construct, condition, seed, source_path.
- `full181_interface_contacts.csv`: the same provenance columns; zero rows
  for a no-contact seed are valid as long as the table/header and seed summary
  exist. Matching rows must refer to the accepted seed-1 source PDB.

Both the summary source and contact sources are checked against the accepted
PDB. A missing analyzer, failing CLI, missing table, duplicate seed summary, or
wrong source blocks seeds 2–15. The root agent should reconcile this CLI contract
with the analysis agent's final entry point before any authorized future run.
No successful contact-chain execution is being claimed in this source review.

## Current execution boundary

Integration update by root: the analyzer is now implemented and accepts
`--condition`, `--seed`, and `--output`. Formal seed-1 execution now passes these
options and writes to its own raw task directory's `analysis/`, avoiding shared
analysis-file races between two GPU branches. Thirteen existing legacy original
PDBs were recalculated with the same generic construct-aware contact/geometry
function; all saved epitope/paratope/CDR counts agreed. This is actual old-data
analysis, not a Full formal validation or a model/smoke run. Full remains 0/30.

The recorded server directory `/home/xinlab/home/ws` is owned by `ws` and has
mode 0700; the current account cannot read the original phase2 environment,
checkpoint, and cache. This remains the reason GPU prediction is unexecuted.
The current release is the baseline. Missing old Native/Immunogen raw ensembles
are not a prerequisite for Full181 prediction and are not requested by this
runner. Existing legacy FASTA/annotations are read for identity checks only.

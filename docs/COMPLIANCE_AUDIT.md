# Attachment 5 compliance audit

Audit date: 2026-09-29. Scope: compact release only. No new GPU prediction was run and no scientific result, sequence, PDB coordinate or experimental value was changed.

## Verification summary

- Quick reconstruction: pass; two runs produced identical SHA-256 `ecaf0fe973cccecc0b43cf51e53641784a8c1200726d22ca830179df01d04326`.
- Strict validator: pass; 10 candidates, 73 indexed logs, six model assets and 26 data files checked.
- Correct environment preflight: pass using the read-only recorded Python 3.10/CUDA environment, verified model cache and HMMER; no output directory was created.
- Default minimal-shell preflight: expected nonzero with explicit missing dependency and model messages; no output directory was created.
- Child failure propagation: pass; synthetic child exit 23 was returned as exit 23 by the full entry.
- Output isolation: pass; all formal steps use `PHASE2_OUTPUT_ROOT`, and preflight with a named output did not create a partial directory.
- Package hashes: pass with `sha256sum -c RELEASE_SHA256SUMS.txt`.
- Absolute path and temporary-artifact scans: pass.
- Full GPU recomputation: not executed in this audit because the instructions prohibit a new GPU run without separate authorization.

## Attachment 5 matrix

| Requirement | Status | Evidence or limitation |
|---|---|---|
| Reproducible runtime and dependencies | Pass | Exact base and full dependency files, CUDA/hardware record and strict version preflight are present. |
| Standard project structure | Pass | README, data, source, models, results, logs, licenses and documentation are present. |
| Single result-generation entry | Pass | `run.sh` deterministically reconstructs the final table. |
| Complete algorithm demonstration | Pass with execution pending | `run_full.sh` invokes the actual formal Boltz workflow and passes preflight in the recorded environment; a new GPU inference run was not authorized. |
| Training entry and training logs | Not applicable | The project did not train or fine-tune the foundation models. |
| Executable Notebook | Not applicable | Attachment 5 treats a Notebook as encouraged; documented quick and full entry points cover the workflows. |
| Data source license and split disclosure | Pass | Every data file has a provenance row; no hidden test set is used. |
| Model versions weights and Model Card | Pass | Six assets have verified size/hash/source records; runtime loading paths are checked. |
| Open-source and external platform disclosure | Pass | Versions, licenses, MSA source and usage limits are documented. |
| Standard final result file | Pass | UTF-8 `results.csv` contains 10 candidates, stable fields, metrics, ranking basis and notes. |
| Structure file conventions | Pass | Coordinate units, chains, numbering, glycosylation, relaxation and model limitations are documented and validated. |
| Portable paths and failure behavior | Pass | Formal steps use isolated relative outputs; preflight and child failures return nonzero. |
| Intellectual property notices | Pass | Boltz and ImmuneBuilder licenses are included; paper PDFs with unverified redistribution permission were removed. |
| Authenticity and reproducibility verification | Pass with execution pending | Quick reconstruction and all static/runtime preflight checks pass; full GPU output regeneration remains unexecuted under the no-GPU instruction. |

Overall conclusion: the release package has no known correctable Attachment 5 defect after this audit. Complete end-to-end GPU result regeneration remains an execution item requiring separate authorization; therefore this audit does not claim that a new full GPU reproduction was completed.

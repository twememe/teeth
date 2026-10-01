# Full181 supplement progress

## 2026-09-21

- Complete source prompt read; exact direct-request exclusions recorded in task_plan.md.
- Three agents dispatched for mapping, interface analysis, and baseline/ensemble comparison. Root handles original server and execution.
- Local compact root: `D:/IL-24(teeth）/phase2_release`; old model directory contains README only. No AGENTS.md found in authoritative project.
- Existing SSH alias `gpu-server` connects to `xinlab-Super-Server`, account pxh. Two A6000 GPUs idle at first audit. Home filesystem has only approximately 770 MiB free; /mnt has approximately 1.6 TiB free. Original project location still being located; no relocation performed.
- Fixed formal schedule: Full181 × {empty, formal_msa} × seeds 1–15; no pilot/smoke model runs, no old GPU tasks.

## Audit and actual MSA completion

- Direct sequence verification passed: upload lines 60/60/60/1; Full/PDB exact 181; H180/L181 preserved; Native155/Immunogen134 are the defined slices; antibody 119/106 and validated IMGT/CDR annotations agree.
- Numbering branch completed 9 requested contract tests, exit 0; no model or hash tests.
- User clarified that Phase2 was run by a collaborator and delivered locally.
- `docs/pip_freeze_final.txt` identifies the historical server project as `/home/xinlab/home/ws/phase2`. On the existing `gpu-server`, `/home/xinlab/home/ws` has mode 0700 and owner/group ws. The current pxh account receives Permission denied for the project, venv, and boltz1 checkpoint. No sudo, account switching, chmod, or permission workaround attempted.
- Formal GPU prediction is blocked at original-resource access, before any seed starts. Empty and formal_msa remain planned 15 each, actual 0 each. The /home filesystem is also near full; /mnt has room but `/mnt/pxh` itself is not writable by this account. No writes to unrelated projects.
- Full A/H/L grouped ColabFold/MMseqs2 MSA completed, 2026-09-21 08:53:17–08:54:35 UTC. Both paired and unpaired raw products, API responses, queries, final CSVs, and process logs are preserved.
- MSA tool is the official Boltz 2.2.1 source downloaded as a source-only wheel; no Boltz model/environment was installed or replaced. Source retrieval initially rejected the workstation Python 3.13 compatibility tag; explicit Python 3.10 wheel targeting resolved this without installing packages.
- Local WSL environment from old Phase1 notes does not exist at its recorded path on this workstation. It was not recreated. Current mapping/analysis preparation uses the existing Windows Python.
- A manifest inspection initially used Windows' default text encoding and failed JSON decoding; explicit UTF-8 reading succeeds. This affected inspection only, not the valid UTF-8 output or MSA search.
- User confirmed the collaborator delivered only the current Phase2 folder and no per-seed raw results. Treat missing legacy originals as a fixed limitation; do not request them again or make them a precondition for new Full tasks. Current GPU block remains original environment/checkpoint access plus suitable writable storage, separately from legacy-results absence.

## Completed independent work and blocked handoff

- Baseline agent resumed after a transport failure and delivered the actual ensemble exporter: 30 saved formal-MSA models, 1899 saved contact rows, 15 PDB files / 13 distinct old models. Saved contact checks 30/30 and saved summary regression 22/22 agree.
- Interface agent suffered a second transport failure without a persisted analyzer. Root implemented the scoped analyzer in two small edits and executed actual recalculation of 13 distinct existing original PDBs. Epitope sets/counts, paratope and CDR counts all agree with saved values; geometry reported only for Native 11 and Immunogen 2 available structures.
- Runner review fixed strict A/H/L sequence validation, finite PDB/PAE/confidence checks, UTF-8 reads, original-version checks, grouped MSA direct CSV reconstruction, formal seed-1 downstream gating and per-task analysis output isolation. No formal GPU/model test executed.
- Current data exports: 30 Full missing-state summaries, 0 Full contact rows, 30 Full QC and PRODIGY status rows, 30 Full rank rows, 6 unselected Full representative slots, three-construct and common-region long tables. No Full PDB/confidence/PAE files were created.
- Final requested numbering run: 9 tests passed. Delivery metadata audit: 30 unique task keys, 0 completed, 30 blocked; 96 old files checked by size/mtime, zero metadata changes. No hash was computed or claimed invariant.
- Acceptance report and plotting handoff are under reports/full181_supplement. Report explicitly states incomplete formal prediction and separate missing legacy/prior limitations.
- Independent source review found a possible half-written shared manifest during two-branch parallel operation. Root corrected this with an OS lock and atomic CSV publication, and atomic JSON state publication. The real blocked-task manifest export succeeds; no parallel GPU validation is claimed.
- Both Full YAML contracts and reconstruction of all three final MSA CSVs from this batch's retained paired/unpaired A3M were actually validated with the final runtime input checker. This was direct input audit, not a new search, model test, or smoke run.
- Final evidence scopes remain separate: 9 numbering tests, 30 saved contact/epitope checks, 22 saved-summary regressions, 13 actual old-PDB subset recalculations, and 0 Full predictions. All 13 old PDB subset observations carry at least one legacy geometry warning; models are retained, and no complete 15-model geometry count is inferred.

## 2026-09-21 root-authorized continuation
- User explicitly supplied root credentials for the same existing server. Root authentication succeeded; no credentials saved in project files.
- Original project /home/xinlab/home/ws/phase2, original venv and Boltz1 checkpoint/CCD are present. Boltz 2.2.1, torch 2.6.0+cu124, PRODIGY 2.4.0, Python 3.10.12 match historical records. Two A6000s idle.
- /home has only 770 MiB free. Created isolated supplement workspace /mnt/il24_full181_supplement_20260921 and original-project link full181_supplement_work on the same server. Original data/results/cache are reused as read-only inputs. Old scripts and outputs are not overwritten.
- Saved baseline size/mtime inventory for 1166 original server files before formal execution. No hash computed.
- Original frozen phase1_prior_deduplicated.csv and old formal empty/MSA raw directories are now readable. Previous local compact-release missing-material conclusions remain historical; new evidence will supersede them after actual validation.
- Existing Full inputs and grouped MSA were transferred, without rerunning searches. Each branch seed1 remains a full 200-step formal task counted among the 30.
- Both formal seed1 tasks completed actual prediction, contact/CDR/QC and PRODIGY before their respective seed2. Two A6000 lanes running concurrently, with fixed formal parameters and no extra model trials.
- User additionally requested consolidated server storage /mnt/pxh/teeth with phase1ab and phase2. Original phase2 was copied there preserving the original; live Full paths remain stable until completion and will be copied into the consolidated phase2. Phase1ab source is the local top-level project excluding phase2_release; WSL tar preserves Linux symlinks in its existing environment/work files. No Phase1 computation launched.
- Source-only baseline adapter updated for both complete server conditions and recovered frozen prior. Agent completed actual validation of all 60 original legacy PDB/confidence/PAE bundles on CPU: all exact sequences, confidence values, explicit chain-pair maps, epitope/paratope/CDR and contact keys agree. Of 4459 contacts, two saved 3-decimal distances differ at the last rounding digit; maximum absolute difference 0.00050029 A, below saved precision tolerance, no threshold changes.
- All 60 old raw models carry legacy atom-clash geometry warnings but no chain breaks, missing backbone or peptide-CN distance outliers; none excluded. New files are in legacy_raw_audit with per-model origins.
- At last observed checkpoint, Full formal accepted 20/30 (empty 11, formal_msa 9), zero failures; both GPUs continued. Interim analysis uses only actual accepted outputs and will be finalized once all 30 finish.

## Final accepted scientific results
- Full 30/30 completed and accepted on server, then 30/30 copied artifacts semantically validated locally; no failed or duplicate model runs.
- Contacts 2232 rows, CDR/FR 420 rows, PRODIGY 30/30 numeric, 90 distinct cross-construct model records, 1808 common-region comparison rows. Final representatives: empty2/6/11 and formal_msa6/9 (seed9 holds two titles).
- Both ensemble exports and actual processed MSA audits pass. 44 legacy summary regressions pass; original server1166/local96 metadata unchanged.
- Native SFTP resumed the preserved phase1ab partial upload after acceptance-package retrieval; archive3075735497 bytes transferred, tar exit0. Source and extracted ordinary file counts24358 and byte total7153573140 agree;2427 links preserved. Original phase2 copy29947 files/4 links agrees by size/mtime/link target.

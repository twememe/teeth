Completed on 2026-09-22: four full LoRA epochs, 24,996 optimizer steps, patience-3 early stopping. Best saved model is epoch 3 / step 14,000. Final common-validation MCC: frozen 0.385864, LoRA 0.499083 (measured on the validation split that also performed checkpoint selection; not an independent confirmation). All saved-model and independent metric checks passed; see output/final_report.md and output/final_independent_audit.json. The historical recovery instructions below are records, not a request to restart.

# Antibody LoRA feasibility experiment

## Capability status (verified 2026-09-30)

- **Implemented and measured**: a single-sequence binary classifier prototype - ESM-2 650M with
  LoRA on `query`/`value` of blocks 0-15 plus a linear head - trained and evaluated on pooled
  public antibody datasets. Model input is exactly one string, `model_sequence`.
- **Not implemented**: antigen-conditioned scoring, wild-type/mutant comparison under one antigen,
  any affinity or mutation-effect prediction, and any feedback loop from IL-24 experimental data.
  There is no code interface that accepts an antigen or a wild-type sequence; the `wildtype`,
  `site` and `mutation` columns that exist in the data tables are provenance for the AbAgym source
  and are never read by the model.
- The reported MCC (frozen 0.385864 -> LoRA 0.499083, delta 0.113219) is measured on the validation
  split that also performed checkpoint selection, so it is not an independent confirmation.
- Details and evidence: `logic_closeout/20260930/closeout.md`.

All work is under this directory (`inovate/`); `scripts/` derives the project root from its own location, and `training_config.yaml` stores paths relative to it, so the workflow does not depend on an absolute repository path.

Start the complete, deadline-aware run with `.venv/bin/python scripts/run_pipeline.py`. The controller takes an exclusive lock to prevent duplicate work. Do not start a second job while `pipeline_status.json` reports an active controller/child PID.

Read `progress.txt`, `pipeline_status.json` and `logs/controller.log` for progress. Failures and repairs are recorded in `error_log.txt`. The absolute deadlines are in `run_state.json`; restarting must never reset the 13-hour budget.

The saved model is an adapter plus linear head, restored against `cache/esm2_650m_hf` using its recorded SHA256 identity. Final model: `output/esm2_650m_lora_best.pt`. Final report: `output/final_report.md`. Those final files are only available after their respective measured phases succeed. Small smoke results use the `smoke_` prefix and are not final performance evidence.

Training defaults: ESM-2 650M; layers 0-15 query/value LoRA; rank 8, alpha 16, dropout 0.1; learning rate 1e-4; global batch 32 (two GPUs, 16 each); at most 15 epochs; patience 3 complete epochs; primary selection metric MCC. The frozen linear-head baseline is fit efficiently on cached embeddings. See YAML for its separately documented optimizer schedule.

Data preparation: `scripts/download_data.py`, then `scripts/preprocess.py`. Raw data versions and hashes are in `data/download_manifest.json`. Exact-sequence groups are isolated between partitions. Biological label direction, mutant partner, assay and representation remain in the CSV metadata. Mixed labels and repeated sequences with assay-dependent outcomes limit biological interpretation.

Dependencies: `.venv` stores all newly installed packages. The pre-existing Python 3.10 CUDA PyTorch/scientific runtime is reused read-only through an explicit `.pth`; it is not modified. New caches and temporary files remain inside this directory. PEFT + Transformers is the user-authorized fallback after the Ab-Tune pinned torch dependency conflict.

## Training recovery on 2026-09-22

At global step 10000, validation created cached rotary tensors under `torch.inference_mode()` and the following same-length training batch raised an autograd error. The failure and original code/checkpoints are preserved under `output/recovery_20260922_0520`. The single targeted repair changes reusable-model evaluation to `torch.no_grad()`; 22 relevant tests and two real 650M two-GPU, 1006-token evaluation/backward transitions passed.

The current configuration continues from the immutable step-10000 checkpoint. Already consumed batches are skipped, prior validation-best weights and completed-epoch patience history are retained, and the original absolute deadlines remain unchanged. The original checkpoint did not save AdamW or RNG state, so these states restart; the final report explicitly discloses this limitation. The repaired production run is the only allowed training retry. Do not launch a duplicate controller or reset its attempt count.

Live status: `pipeline_status.json`; repair evidence: `output/training_repair_report.json`, `output/repair_gpu_probe.json`, `logs/training_repair_tests.log`; current training log: `logs/lora_training_attempt1.log`.

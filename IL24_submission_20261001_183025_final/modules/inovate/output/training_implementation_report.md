# Training implementation handoff

Implemented scripts/train.py, training_config.yaml, and tests/test_training.py. No real GPU training was launched by the implementation worker. Ten CPU tests pass; syntax compilation and CLI help pass. The tests verify complete ordered extraction stitching, rejection of partial caches, exact validation sharding, equal DDP optimization batches, strict sequence-length checks, metrics including constant-label behavior, patience, missing checkpoint tensor rejection, and a tiny synthetic one-epoch adapter/head training loop with checkpoint and prediction exports. Temporary test files stay under inovate/tmp.

## Commands from /mnt/pxh/teeth/inovate

```
.venv/bin/python -m torch.distributed.run --standalone --nproc_per_node=2 scripts/train.py --mode smoke --config training_config.yaml --max-seconds 2400
.venv/bin/python -m torch.distributed.run --standalone --nproc_per_node=2 scripts/train.py --mode extract --config training_config.yaml --max-seconds 7200
.venv/bin/python scripts/train.py --mode baseline --config training_config.yaml --max-seconds 1800
.venv/bin/python -m torch.distributed.run --standalone --nproc_per_node=2 scripts/train.py --mode lora --config training_config.yaml --max-seconds 32400 --deadline 2026-09-22T12:37:56.557353+08:00
.venv/bin/python -m torch.distributed.run --standalone --nproc_per_node=2 scripts/train.py --mode evaluate --config training_config.yaml --max-seconds 2400
```

Pass the absolute experiment cutoff with --deadline to every phase if orchestrator time is tighter than its suggested phase limit. --max-seconds and --deadline use the earlier limit. Test data remains sealed unless an explicit --mode evaluate --split test command is issued. No retry loop exists.

## Implemented experiment

ESM-2 650M model construction belongs to scripts/modeling.py. Initial config retains gradient checkpointing true and SDPA. Adapter layers 0 through 15, query/value projections, r=8, alpha=16, dropout=0.1. AdamW lr=1e-4, weight decay 0.01, 15 epochs maximum, patience 3 completed epochs on validation MCC, gradient clipping 1.0, BF16 CUDA autocast with FP32 model parameters. DDP requires two ranks for formal training; per-rank batch 16 gives global batch 32. Model sequences are actual model_sequence CSV values. Tokenization forbids silent truncation and verifies one residue per token.

Frozen baseline uses distributed frozen feature extraction with disjoint sample indices, then one identical Linear(1280,1) head fitted efficiently on rank 0. Baseline optimizer deliberately differs: AdamW lr=1e-3, batch 4096, up to 30 epochs, patience 3. LoRA warm starts from that fitted frozen head; this is explicitly saved in checkpoint training metadata. Source checkpoint SHA256 identities, baseline checkpoint hash, precise config, epoch, global steps, and metrics accompany selected weights.

Long sequences are length bucketed. Training drops the incomplete global batch (31 rows with current 199999-row train split) for equal DDP steps; validation and extraction have no duplicate padding and stitch exact original indices. Deadlines are collectively checked before training batches. Extraction/evaluation ranks can stop independently and then gather status at the end, avoiding per-batch collectives across unequal shard lengths. Incomplete extraction and incomplete standalone evaluation exit with errors. If the training deadline arrives after optimizer steps but without a complete validation-selected checkpoint, lora_run.json reports incomplete_needs_final_validation and the latest_checkpoint path; the controller can evaluate that checkpoint in its reserved final phase. Zero optimizer steps raise an error. Last adapter/head state is saved separately at checkpoints and training stop. Optimizer state is not saved and automatic resumption is not implemented.

Validation selection occurs every 2000 optimizer steps and after each epoch. Best selection can occur at an interval; early-stopping patience counts completed epochs. Evaluate only complete validation passes. At a deadline, the existing best remains usable and the latest trainable weights are retained; a partial validation pass never produces metrics. Actual started/completed epochs and steps are written to lora_run.json.

Measured random baseline uses seeded independent uniform probabilities and the same validation rows. Theoretical balanced F1=0.5 and MCC=0 remain a separate reference. MCC, F1, and Spearman(probability, binary label) are reported overall, by source, and by source/mutated partner. Within-source raw-score Spearman is additional descriptive output; different assay directions may differ. Labels are described only as numerical high-score labels, never as a shared biological positive meaning.

## Outputs and fallback controls

Formal: output/esm2_650m_lora_best.pt, esm2_650m_lora_latest.pt, esm2_650m_frozen_best.pt, losses.csv, validation_predictions.csv, per-model prediction CSV and metrics JSON, metrics.json, comparison.csv, comparison.md, baseline_run.json, lora_run.json. Smoke weights, caches, losses, metrics and comparisons all use smoke_ prefixes. Smoke metrics are kept out of formal comparisons.

Progress and full error logs are appended both under output/ and at ROOT/progress.txt, ROOT/error_log.txt; per-rank logs live under output/. Root orchestrator should redirect torchrun stdout/stderr to a stage log as usual. Config supports gradient_accumulation_steps (default 1); an explicitly authorized OOM fallback to per_device_batch_size=8 plus gradient_accumulation_steps=2 preserves global batch 32. extraction_batch_size and evaluation_batch_size can be reduced separately. Trainable checkpoints reject missing head or adapter tensors during evaluation. All phases read environment.json cache settings.

## Follow-up review fixes

Embedding reuse compares both the data hash and full source model identity (model directory, configuration SHA256 and weight SHA256). Stale caches are recomputed during extraction and rejected during baseline fitting. LoRA warm start verifies source model identity and both split data identities. Baseline checkpoints written before this schema require refitting the cached linear head before LoRA warm start. Loss rows now include a globally averaged training loss at every log_every_steps interval, in addition to validation/epoch rows. Ten focused CPU tests cover the new deadline handoff, zero-step failure, model-cache identity rejection, warm-start identity rejection, and loss-history cadence. No GPU work was launched for these changes.

#!/usr/bin/env python
"""Run official standalone BepiPred-3.0 vt_pred and standardize its output."""

from __future__ import annotations

import argparse
import csv
import math
import os
import statistics
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


DEFAULT_THRESHOLD = 0.1512
ROLLING_WINDOW = 9


def read_single_fasta(path: Path) -> tuple[str, str]:
    header = ""
    sequence_parts: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                raise ValueError(f"Expected one FASTA record in {path}")
            header = line[1:]
        else:
            sequence_parts.append(line)
    if not header or not sequence_parts:
        raise ValueError(f"Invalid or empty FASTA: {path}")
    return header, "".join(sequence_parts)


def _run_command(command: list[str], env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=env,
        check=False,
    )


def build_official_command(
    bp3_python: Path,
    cli: Path,
    input_fasta: Path,
    raw_dir: Path,
) -> list[str]:
    return [
        str(bp3_python),
        str(cli),
        "-i",
        str(input_fasta),
        "-o",
        str(raw_dir),
        "-pred",
        "vt_pred",
        "-esm_dir",
        str(raw_dir / "esm_encodings"),
    ]


def run_official_cli(
    bp3_python: Path,
    repo: Path,
    input_fasta: Path,
    raw_dir: Path,
    log_path: Path,
) -> tuple[str, float]:
    cli = repo / "bepipred3_CLI.py"
    if not cli.is_file():
        raise FileNotFoundError(f"Official CLI not found: {cli}")
    raw_dir.mkdir(parents=True, exist_ok=True)
    command = build_official_command(bp3_python, cli, input_fasta, raw_dir)
    base_env = os.environ.copy()
    attempts: list[tuple[str, subprocess.CompletedProcess[str]]] = []

    started = time.perf_counter()
    first = _run_command(command, base_env)
    attempts.append(("CUDA-preferred official default", first))
    device = "cuda"

    oom_text = first.stdout.lower()
    if first.returncode != 0 and (
        "cuda out of memory" in oom_text or "cuda error: out of memory" in oom_text
    ):
        cpu_env = base_env.copy()
        cpu_env["CUDA_VISIBLE_DEVICES"] = ""
        second = _run_command(command, cpu_env)
        attempts.append(("CPU retry after explicit CUDA OOM", second))
        first = second
        device = "cpu_after_cuda_oom"

    inference_seconds = time.perf_counter() - started

    log_parts = [
        f"utc_started={datetime.now(timezone.utc).isoformat()}",
        "official_repo=https://github.com/UberClifford/BepiPred-3.0",
        f"command={' '.join(command)}",
        "model_mode=vt_pred",
        f"official_default_threshold={DEFAULT_THRESHOLD}",
        f"official_default_rolling_window={ROLLING_WINDOW}",
        f"runtime_device={device}",
        f"inference_seconds={inference_seconds:.3f}",
    ]
    for index, (label, completed) in enumerate(attempts, start=1):
        log_parts.extend(
            [
                "",
                f"=== ATTEMPT {index}: {label} ===",
                f"returncode={completed.returncode}",
                completed.stdout,
            ]
        )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("\n".join(log_parts), encoding="utf-8")
    if first.returncode != 0:
        raise RuntimeError(
            f"Official BepiPred-3.0 exited {first.returncode}; complete log: {log_path}"
        )
    return device, inference_seconds


def parse_official_outputs(raw_dir: Path, reference_fasta: Path) -> list[dict]:
    _, reference = read_single_fasta(reference_fasta)
    raw_csv = raw_dir / "raw_output.csv"
    _, prediction_sequence = read_single_fasta(raw_dir / "Bcell_epitope_preds.fasta")
    if len(reference) != 181:
        raise ValueError(f"Frozen reference must be 181 aa; found {len(reference)}")

    with raw_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        expected = [
            "Accession",
            "Residue",
            "BepiPred-3.0 score",
            "BepiPred-3.0 linear epitope score",
        ]
        if reader.fieldnames != expected:
            raise ValueError(f"Unexpected official raw_output.csv columns: {reader.fieldnames}")
        official_rows = list(reader)

    if len(official_rows) != 181 or len(prediction_sequence) != 181:
        raise ValueError(
            "BepiPred output must contain exactly 181 residues; "
            f"raw={len(official_rows)}, prediction={len(prediction_sequence)}"
        )

    rows: list[dict] = []
    for position, (raw, project_aa, prediction_aa) in enumerate(
        zip(official_rows, reference, prediction_sequence), start=1
    ):
        raw_aa = raw["Residue"].strip().upper()
        if raw_aa != project_aa or prediction_aa.upper() != project_aa:
            raise ValueError(
                f"Coordinate mismatch at project position {position}: "
                f"reference={project_aa}, raw={raw_aa}, prediction={prediction_aa}"
            )
        probability = float(raw["BepiPred-3.0 score"].strip())
        rolling = float(raw["BepiPred-3.0 linear epitope score"].strip())
        if not math.isfinite(probability) or not math.isfinite(rolling):
            raise ValueError(f"Non-finite BepiPred score at position {position}")
        official_prediction = int(prediction_aa.isupper())
        computed_prediction = int(probability >= DEFAULT_THRESHOLD)
        if official_prediction != computed_prediction:
            raise ValueError(
                f"Official FASTA/default-threshold disagreement at position {position}"
            )
        rows.append(
            {
                "position": position,
                "AA": project_aa,
                "bepipred_positive_probability": probability,
                "bepipred_default_threshold_prediction": official_prediction,
                "bepipred_rolling_mean_9aa": rolling,
                "official_default_threshold": DEFAULT_THRESHOLD,
                "model_mode": "vt_pred",
            }
        )
    return rows


def write_standard_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_profile(path: Path, rows: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    positions = [row["position"] for row in rows]
    probabilities = [row["bepipred_positive_probability"] for row in rows]
    rolling = [row["bepipred_rolling_mean_9aa"] for row in rows]
    fig, axis = plt.subplots(figsize=(11, 4.8))
    axis.plot(positions, probabilities, color="#2f6f9f", linewidth=1.1, label="Positive probability")
    axis.plot(positions, rolling, color="#d95f02", linewidth=1.5, label="Official 9-aa rolling mean")
    axis.axhline(
        DEFAULT_THRESHOLD,
        color="#4d4d4d",
        linestyle="--",
        linewidth=1,
        label=f"Official default threshold ({DEFAULT_THRESHOLD})",
    )
    axis.set(xlabel="Project position (1–181)", ylabel="BepiPred-3.0 score", xlim=(1, 181))
    axis.grid(alpha=0.2)
    axis.legend(frameon=False, ncol=3, fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def summarize_high_regions(rows: list[dict]) -> list[dict]:
    regions = []
    start = None
    probabilities: list[float] = []
    for row in rows + [{"position": len(rows) + 1, "bepipred_default_threshold_prediction": 0}]:
        if int(row["bepipred_default_threshold_prediction"]) == 1:
            if start is None:
                start = int(row["position"])
                probabilities = []
            probabilities.append(float(row["bepipred_positive_probability"]))
        elif start is not None:
            end = int(row["position"]) - 1
            regions.append(
                {
                    "start": start,
                    "end": end,
                    "length": end - start + 1,
                    "mean_probability": statistics.fmean(probabilities),
                }
            )
            start = None
            probabilities = []
    return regions


def write_report(
    path: Path,
    rows: list[dict],
    device: str,
    inference_seconds: float,
    log_path: Path,
) -> None:
    probabilities = [row["bepipred_positive_probability"] for row in rows]
    predicted = sum(row["bepipred_default_threshold_prediction"] for row in rows)
    top = max(rows, key=lambda row: row["bepipred_positive_probability"])
    top_residues = sorted(
        rows, key=lambda row: row["bepipred_positive_probability"], reverse=True
    )[:10]
    high_regions = sorted(
        summarize_high_regions(rows),
        key=lambda region: (-region["mean_probability"], region["start"]),
    )
    top_residue_text = ", ".join(
        f"{row['position']}{row['AA']} ({row['bepipred_positive_probability']:.4f})"
        for row in top_residues
    )
    high_region_text = "; ".join(
        f"{region['start']}–{region['end']} (n={region['length']}, mean={region['mean_probability']:.4f})"
        for region in high_regions[:10]
    ) or "None under the official default threshold"
    report = f"""# Step 3 — BepiPred-3.0

## Step 3 完成

### 1. 本步做了什么
对冻结的 181-aa 项目参考运行官方 standalone BepiPred-3.0，并保存逐 residue 正类概率、官方默认阈值分类和官方 9-aa rolling mean。

### 2. 使用的算法
BepiPred-3.0 official standalone，`vt_pred`，默认阈值 {DEFAULT_THRESHOLD}，默认 rolling window {ROLLING_WINDOW}。

### 3. 算法属于什么
基于预训练 ESM-2 protein language model 表征的深度学习推理；没有训练或微调模型。

### 4. 算法原理
通俗地说，ESM-2 从大量蛋白序列中学习氨基酸“语言规律”，BepiPred 再判断每个位置多像已知抗体接触 residue。技术上，官方五折前馈网络 ensemble 对 ESM-2 residue embedding 输出正类概率，`vt_pred` 在平均正类概率上应用固定阈值。

### 5. 为什么本项目需要它
它提供独立的 1D sequence-based epitope evidence；后续融合应优先使用连续概率，而不是只使用 0/1 分类。

### 6. 输入
- `inputs/reference_181.fasta`

### 7. 输出
- `results/03_bepipred/raw/`
- `results/03_bepipred/bepipred_residue_scores.csv`
- `results/03_bepipred/bepipred_profile.png`
- `reports/step03_bepipred.md`

### 8. 关键数值结果
- Residue rows：{len(rows)}。
- Positive probability：min={min(probabilities):.6f}，mean={statistics.fmean(probabilities):.6f}，max={max(probabilities):.6f}。
- 官方默认阈值阳性 residue：{predicted}/{len(rows)}。
- 最高连续概率位置：project {top['position']} ({top['AA']})，score={top['bepipred_positive_probability']:.6f}。
- 高分 residue（Top 10）：{top_residue_text}。
- 官方阈值连续阳性区域：{high_region_text}。
- Runtime device result：{device}。
- Inference time：{inference_seconds:.2f} seconds。

### 9. 当前结果怎样解释
高分位置是 sequence-based computational evidence，不等于实验验证表位、受体界面或中和位点。9-aa rolling mean 仅用于线性连续性辅助展示。

### 10. 是否存在问题
模型默认阈值不是对本蛋白重新校准的阈值，不能用历史结果调参。完整官方运行日志位于 `{log_path.as_posix()}`。

### 11. 下一步
Integration 将连续 BepiPred 概率与当前可用的 RSA 表面证据和 conservation 证据透明融合。DiscoTope 分支标记为缺失，不引入替代模型。
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(report, encoding="utf-8")


def write_failure_report(path: Path, message: str, log_path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# Step 3 — BepiPred-3.0\n\n"
        "## 状态\n\nOfficial standalone inference failed; BepiPred evidence is missing and was not replaced.\n\n"
        f"## Error\n\n`{message}`\n\n"
        f"Complete log: `{log_path.as_posix()}`\n",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bp3-python", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--scores-csv", type=Path, required=True)
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        device, inference_seconds = run_official_cli(
            args.bp3_python, args.repo, args.input, args.raw_dir, args.log
        )
        rows = parse_official_outputs(args.raw_dir, args.input)
        write_standard_csv(args.scores_csv, rows)
        write_profile(args.profile, rows)
        write_report(args.report, rows, device, inference_seconds, args.log)
    except Exception as error:
        write_failure_report(args.report, str(error), args.log)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

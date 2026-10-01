#!/usr/bin/env python3
"""FreeSASA surface accessibility and DSSP annotation for IL-24 residues 1-181."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import tempfile
import traceback
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from statistics import mean, median

from Bio import __version__ as BIOPYTHON_VERSION
from Bio.PDB import DSSP, PDBParser

try:
    from .structure_qc import analyze_structure, load_reference
except ImportError:  # Direct execution as ``python scripts/surface_analysis.py``.
    from structure_qc import analyze_structure, load_reference


class ToolUnavailableError(RuntimeError):
    """Raised when a required external structural tool cannot be used."""


def _require_complete_mapping(rows: list[dict], label: str) -> None:
    positions = [row["position"] for row in rows]
    if len(rows) != 181 or positions != list(range(1, 182)):
        raise ValueError(f"{label} did not map exactly to project positions 1-181")


def _freesasa_version(module) -> str:
    for attribute in ("__version__", "version"):
        value = getattr(module, attribute, None)
        if callable(value):
            value = value()
        if value:
            return str(value)
    try:
        from importlib.metadata import version

        return version("freesasa")
    except Exception:
        pass
    return "not_exposed_by_python_binding"


@contextmanager
def _freesasa_compatible_path(pdb_path: str | Path):
    """Stage an exact copy when the FreeSASA binding cannot encode the path."""
    source = Path(pdb_path)
    try:
        str(source).encode("ascii")
        yield str(source)
    except UnicodeEncodeError:
        with tempfile.TemporaryDirectory(prefix="il24_freesasa_") as directory:
            staged = Path(directory) / "input.pdb"
            shutil.copyfile(source, staged)
            yield str(staged)


def calculate_freesasa(pdb_path: str | Path, structure_rows: list[dict]) -> tuple[list[dict], dict]:
    """Calculate absolute SASA and FreeSASA-normalized RSA for every residue."""
    try:
        import freesasa
    except ImportError as error:
        raise ToolUnavailableError("FreeSASA Python bindings are not installed") from error

    parameters = freesasa.Parameters()
    with _freesasa_compatible_path(pdb_path) as compatible_path:
        freesasa_structure = freesasa.Structure(compatible_path)
        result = freesasa.calc(freesasa_structure, parameters)
    residue_areas = result.residueAreas()
    output: list[dict] = []

    for structure_row in structure_rows:
        chain_id = structure_row["pdb_chain"]
        residue_number = str(structure_row["pdb_resseq"])
        icode = structure_row["pdb_icode"]
        chain_areas = residue_areas.get(chain_id, {})
        area = chain_areas.get(residue_number)
        if area is None and icode:
            area = chain_areas.get(f"{residue_number}{icode}")
        if area is None:
            raise ValueError(
                f"FreeSASA result is missing chain {chain_id} residue {residue_number}{icode}"
            )
        absolute_sasa = float(area.total)
        rsa = float(area.relativeTotal)
        if not math.isfinite(absolute_sasa) or not math.isfinite(rsa):
            raise ValueError(
                f"FreeSASA returned a non-finite value for project position {structure_row['position']}"
            )
        output.append(
            {
                "position": structure_row["position"],
                "aa": structure_row["aa"],
                "pdb_chain": chain_id,
                "pdb_resseq": structure_row["pdb_resseq"],
                "absolute_sasa_angstrom2": round(absolute_sasa, 6),
                "rsa": round(rsa, 6),
            }
        )

    _require_complete_mapping(output, "FreeSASA")
    algorithm = str(parameters.algorithm())
    metadata = {
        "freesasa_available": True,
        "freesasa_version": _freesasa_version(freesasa),
        "freesasa_algorithm": algorithm,
        "probe_radius_angstrom": float(parameters.probeRadius()),
        "total_sasa_angstrom2": round(float(result.totalArea()), 6),
    }
    return output, metadata


def _resolve_executable(executable: str | Path) -> str:
    candidate = Path(executable)
    if candidate.is_file():
        return str(candidate.resolve())
    resolved = shutil.which(str(executable))
    if resolved is None:
        raise ToolUnavailableError(f"DSSP executable not found: {executable}")
    return resolved


def _secondary_structure_class(code: str) -> str:
    if code in {"H", "G", "I"}:
        return "helix"
    if code in {"E", "B"}:
        return "strand"
    return "coil"


def calculate_dssp(
    pdb_path: str | Path,
    structure_rows: list[dict],
    dssp_executable: str | Path = "mkdssp",
) -> tuple[list[dict], dict]:
    """Run mkdssp and map secondary-structure annotations to project positions."""
    executable = _resolve_executable(dssp_executable)
    parser = PDBParser(QUIET=True)
    model = parser.get_structure("IL24_DSSP", str(pdb_path))[0]
    try:
        dssp = DSSP(model, str(pdb_path), dssp=executable)
    except Exception as error:
        raise ToolUnavailableError(f"DSSP execution failed: {error}") from error

    output: list[dict] = []
    for structure_row in structure_rows:
        residue_id = (
            " ",
            int(structure_row["pdb_resseq"]),
            structure_row["pdb_icode"] or " ",
        )
        key = (structure_row["pdb_chain"], residue_id)
        if key not in dssp:
            raise ValueError(
                f"DSSP result is missing project position {structure_row['position']} ({key})"
            )
        record = dssp[key]
        raw_code = record[2] if record[2] != "-" else "C"
        output.append(
            {
                "position": structure_row["position"],
                "aa": structure_row["aa"],
                "pdb_chain": structure_row["pdb_chain"],
                "pdb_resseq": structure_row["pdb_resseq"],
                "dssp_code": raw_code,
                "secondary_structure_class": _secondary_structure_class(raw_code),
                "dssp_available": True,
            }
        )

    _require_complete_mapping(output, "DSSP")
    version_process = subprocess.run(
        [executable, "--version"],
        capture_output=True,
        text=True,
        timeout=20,
        check=False,
    )
    version_text = (version_process.stdout or version_process.stderr).strip().splitlines()
    metadata = {
        "dssp_available": True,
        "dssp_executable": executable,
        "dssp_version": version_text[0] if version_text else "unknown",
        "biopython_version": BIOPYTHON_VERSION,
    }
    return output, metadata


def write_surface_outputs(
    sasa_rows: list[dict],
    dssp_rows: list[dict],
    metadata: dict,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Write strict per-residue FreeSASA/DSSP tables and the RSA profile."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _require_complete_mapping(sasa_rows, "FreeSASA output")
    _require_complete_mapping(dssp_rows, "DSSP output")
    if [row["aa"] for row in sasa_rows] != [row["aa"] for row in dssp_rows]:
        raise ValueError("FreeSASA and DSSP amino-acid mappings disagree")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sasa_path = output_dir / "residue_sasa.csv"
    dssp_path = output_dir / "residue_secondary_structure.csv"
    plot_path = output_dir / "rsa_profile.png"

    with sasa_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sasa_rows[0]))
        writer.writeheader()
        writer.writerows(sasa_rows)
    with dssp_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dssp_rows[0]))
        writer.writeheader()
        writer.writerows(dssp_rows)

    fig, axis = plt.subplots(figsize=(10, 4.2))
    axis.plot(
        [row["position"] for row in sasa_rows],
        [row["rsa"] for row in sasa_rows],
        color="#26734D",
        linewidth=1.5,
    )
    axis.set(xlabel="Project residue position", ylabel="Relative solvent accessibility (RSA)", xlim=(1, 181))
    axis.set_ylim(bottom=0)
    axis.set_title("IL-24 FreeSASA relative solvent accessibility")
    fig.tight_layout()
    fig.savefig(plot_path, dpi=200)
    plt.close(fig)
    return {
        "residue_sasa": sasa_path,
        "residue_secondary_structure": dssp_path,
        "rsa_profile": plot_path,
    }


def write_surface_report(
    sasa_rows: list[dict], dssp_rows: list[dict], metadata: dict, report_path: str | Path
) -> Path:
    """Write the Step 4 method report from real FreeSASA and DSSP results."""
    sasa = [row["absolute_sasa_angstrom2"] for row in sasa_rows]
    rsa = [row["rsa"] for row in sasa_rows]
    classes = Counter(row["secondary_structure_class"] for row in dssp_rows)
    text = f"""# Step 4：Structural Surface Accessibility

## 1. 本步做了什么

在通过 181-aa 序列硬门控的 AlphaFold 坐标上，逐 residue 计算 absolute SASA 和 RSA，并用 DSSP 添加二级结构注释。所有表严格使用 project 1–181 编号。

## 2. 使用的算法

FreeSASA {metadata['freesasa_version']}（{metadata['freesasa_algorithm']}，probe radius {metadata['probe_radius_angstrom']:.2f} Å）；{metadata['dssp_version']}，通过 Biopython DSSP {metadata['biopython_version']} 调用。

## 3. 算法属于什么

FreeSASA 是几何表面积计算；DSSP 是传统结构生物信息学注释。两者都不是 AI 模型。

## 4. 算法原理

FreeSASA 模拟一个水分子尺度的探针沿蛋白表面滚动，计算探针可接触的原子表面积；absolute SASA 是面积，RSA 是按 residue 类型最大可接触面积归一化后的相对值。DSSP 根据主链氢键与几何模式判定 helix、strand 和 coil。

## 5. 为什么本项目需要它

抗体必须物理接触抗原表面，因此 RSA 提供独立的结构可接近性证据。DSSP 只帮助解释局部结构背景，不自动加分或减分。

## 6. 输入

`inputs/AF-Q925S4-F1-model_v6.pdb`；`inputs/reference_181.fasta`，且结构序列必须完全一致后才运行。

## 7. 输出

`results/04_surface/residue_sasa.csv`、`residue_secondary_structure.csv`、`rsa_profile.png`。

## 8. 关键数值结果

- FreeSASA total SASA：{metadata['total_sasa_angstrom2']:.2f} Å²。
- residue absolute SASA minimum/mean/median/maximum：{min(sasa):.2f} / {mean(sasa):.2f} / {median(sasa):.2f} / {max(sasa):.2f} Å²。
- RSA minimum/mean/median/maximum：{min(rsa):.4f} / {mean(rsa):.4f} / {median(rsa):.4f} / {max(rsa):.4f}。
- DSSP helix/strand/coil：{classes.get('helix', 0)} / {classes.get('strand', 0)} / {classes.get('coil', 0)} residues。

## 9. 当前结果怎样解释

较高 RSA 表示 AlphaFold 模型中更暴露、更容易被溶剂探针接触；它支持表面可接近性，但不能单独证明抗体结合。DSSP 仅为 annotation，不进入最终 score。

## 10. 是否存在问题

SASA/RSA 依赖当前单体 AlphaFold 构象，并受低 pLDDT 区域坐标不确定性影响；低 pLDDT residue 的结构指标必须降级解释。

## 11. 下一步

Integration Agent 仅在严格 181 行映射下使用 RSA；DSSP 保留为解释性注释。
"""
    report_path = Path(report_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(text, encoding="utf-8")
    return report_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--dssp", default="mkdssp")
    parser.add_argument("--failure-log", type=Path, required=True)
    args = parser.parse_args()

    reference = load_reference(args.reference)
    structure_rows, qc = analyze_structure(args.pdb, reference)
    if not qc["analysis_allowed"]:
        raise SystemExit("Structure sequence/completeness gate did not pass")

    try:
        sasa_rows, freesasa_metadata = calculate_freesasa(args.pdb, structure_rows)
        dssp_rows, dssp_metadata = calculate_dssp(args.pdb, structure_rows, args.dssp)
    except Exception:
        args.failure_log.parent.mkdir(parents=True, exist_ok=True)
        args.failure_log.write_text(traceback.format_exc(), encoding="utf-8")
        raise

    metadata = {**freesasa_metadata, **dssp_metadata}
    write_surface_outputs(sasa_rows, dssp_rows, metadata, args.output_dir)
    write_surface_report(sasa_rows, dssp_rows, metadata, args.report)
    args.failure_log.parent.mkdir(parents=True, exist_ok=True)
    args.failure_log.write_text(
        json.dumps({"status": "success", **metadata}, indent=2) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""V2 final within-band structural veto and maximal-retention decision."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


MODULE_ROOT = Path(__file__).resolve().parents[1]
BANDS = MODULE_ROOT / "results/02_hydrophobic_boundary/hydrophobic_boundary_bands.csv"
RISK = MODULE_ROOT / "results/03_structural_boundary/structural_risk.csv"
NTERM = MODULE_ROOT / "results/01_n_terminal/n_terminal_boundary.csv"
SENSITIVITY = MODULE_ROOT / "results/03_structural_boundary/contact_quantile_sensitivity.csv"
HYDRO_ENDPOINTS = MODULE_ROOT / "results/02_hydrophobic_boundary/hydropathy_endpoint_complete_windows.csv"


def select_final_boundary(
    bands: pd.DataFrame,
    risk: pd.DataFrame,
    construct_start: int,
) -> tuple[dict[str, object], pd.DataFrame]:
    if bands.empty:
        return {
            "status": "NO_STABLE_C_TERMINAL_BOUNDARY",
            "start": construct_start,
            "recommended_end": None,
            "recommended_construct": None,
        }, pd.DataFrame()
    if "hydropathy_rank" in bands:
        primary_rows = bands.loc[bands["hydropathy_rank"] == 1]
    else:
        primary_rows = bands.loc[bands["is_primary_band"].astype(bool)]
    if len(primary_rows) != 1:
        raise ValueError("expected exactly one primary hydropathy band")
    primary = primary_rows.iloc[0]
    band_start = int(primary["band_start"])
    band_end = int(primary["band_end"])
    members = pd.DataFrame({"endpoint": list(range(band_start, band_end + 1))})
    members = members.merge(risk, on="endpoint", how="left", validate="one_to_one")
    if members["structural_risk"].isna().any():
        raise ValueError("missing structural risk for a primary-band endpoint")
    members.insert(0, "band_id", primary["band_id"])
    members.insert(1, "hydropathy_rank", int(primary.get("hydropathy_rank", 1)))
    members["is_structurally_acceptable"] = (
        members["structural_risk"] != "HIGH_STRUCTURAL_RISK"
    )
    acceptable = members.loc[members["is_structurally_acceptable"], "endpoint"].astype(int).tolist()
    if not acceptable:
        decision = {
            "status": "STRUCTURAL_CONFLICT",
            "start": construct_start,
            "primary_hydrophobic_boundary_band": [band_start, band_end],
            "acceptable_endpoints": [],
            "recommended_end": None,
            "recommended_construct": None,
        }
        members["is_recommended"] = False
        return decision, members
    recommended_end = max(acceptable)
    members["is_recommended"] = members["endpoint"] == recommended_end
    decision = {
        "status": "RECOMMENDED",
        "start": construct_start,
        "primary_hydrophobic_boundary_band": [band_start, band_end],
        "primary_peak_endpoint": int(primary["peak_endpoint"]),
        "primary_peak_consensus": float(primary["peak_consensus"]),
        "primary_band_mean_consensus": float(primary["mean_consensus"]),
        "acceptable_endpoints": acceptable,
        "recommended_end": recommended_end,
        "recommended_construct": f"{construct_start}-{recommended_end}",
        "construct_length_aa": recommended_end - construct_start + 1,
        "decision_rule": "rightmost structurally acceptable endpoint inside the primary hydropathy band",
        "length_role": "final within-band tie-break only",
        "legacy_global_selector_used": False,
        "phase1b_modified": False,
    }
    return decision, members


def _save_figure(table: pd.DataFrame, decision: dict[str, object], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    colors = table["structural_risk"].map(
        {"ACCEPTABLE_STRUCTURAL_RISK": "#2a9d8f", "HIGH_STRUCTURAL_RISK": "#d00000"}
    )
    ax.scatter(table["endpoint"], [1] * len(table), c=colors, s=120, edgecolor="black")
    for row in table.itertuples():
        ax.text(row.endpoint, 1.035, str(row.endpoint), ha="center", fontsize=9)
    if decision.get("recommended_end") is not None:
        endpoint = int(decision["recommended_end"])
        ax.scatter([endpoint], [1], marker="*", s=420, color="#f9c74f", edgecolor="black", label="selected endpoint")
    ax.set_ylim(0.8, 1.22)
    ax.set_yticks([])
    ax.set_xlabel("Endpoint in primary hydropathy band")
    ax.set_title("Hydrophobic stable band → structural veto → maximal-retention endpoint")
    ax.text(0.01, 0.1, "green = acceptable; red = high structural risk", transform=ax.transAxes)
    ax.legend(frameon=False, loc="lower right")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def _markdown_table(frame: pd.DataFrame) -> str:
    values = frame.astype(str)
    header = "| " + " | ".join(values.columns) + " |"
    separator = "| " + " | ".join(["---"] * len(values.columns)) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in values.to_numpy().tolist()]
    return "\n".join([header, separator, *rows])


def _write_report(decision: dict[str, object], table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    acceptable = ", ".join(map(str, decision.get("acceptable_endpoints", []))) or "none"
    nterm = pd.read_csv(NTERM).iloc[0]
    sensitivity = pd.read_csv(SENSITIVITY)
    structural = pd.read_csv(RISK)
    primary_structural = structural.loc[structural["endpoint"].between(154, 160), [
        "endpoint", "local_cross_contact_density", "dssp_status", "boundary_mean_RSA",
        "boundary_median_pLDDT", "structural_risk",
    ]].copy()
    sensitivity_view = sensitivity[["endpoint", "q60", "q65", "q70", "q75", "q80", "q85", "q90", "acceptable_fraction", "sensitivity_level"]].copy()
    sensitivity_view = sensitivity_view.replace({"ACCEPTABLE_STRUCTURAL_RISK": "OK", "HIGH_STRUCTURAL_RISK": "HIGH"})
    hydro = pd.read_csv(HYDRO_ENDPOINTS)
    peak = hydro.loc[hydro["endpoint"] == int(decision["primary_peak_endpoint"])].iloc[0]
    scale_rows = pd.DataFrame({
        "window": [7, 9, 11, 15, 21],
        "DeltaH_at_peak_158": [f"{float(peak[f'DeltaH_{w}']):.4f}" for w in (7, 9, 11, 15, 21)],
        "percentile": [f"{float(peak[f'P{w}']):.4f}" for w in (7, 9, 11, 15, 21)],
    })
    path.write_text(
        "# IL-24 High-level Phase 1A v2.1\n## Recombinant Immunogen Construct Selection\n\n"
        "## 1. Executive Summary\n\n"
        f"Recommended construct: **{decision.get('recommended_construct')}** ({decision.get('construct_length_aa')} aa). Start 27 is supported by reviewed UniProt, DeepSig SP 1-26 (confidence 1.0), and N-terminal hydropathy. End 160 is the rightmost q75-acceptable endpoint in the stable hydrophobic band 154-160. Endpoint 160 sensitivity is **{decision.get('structural_sensitivity')}** ({decision.get('endpoint_160_acceptable_fraction'):.3f} acceptable).\n\n"
        "## 2. Input integrity\n\nFASTA and PDB Chain A are 181 residues with exact 181/181 identity and project numbering 1-181.\n\n"
        "## 3. N-terminal signal peptide\n\n"
        "- SignalP 6: not executed; `UNAVAILABLE_EXTERNAL_LICENSED_PACKAGE` (official portable package and weights absent).\n"
        f"- Reviewed UniProt Q925S4: project signal peptide 1-26; cleavage {nterm['cleavage']}.\n"
        f"- DeepSig 0.9: euk mode; SP {nterm['deepsig_signal']}; cleavage {nterm['deepsig_cleavage']}; confidence {nterm['deepsig_confidence']}; mature Chain {nterm['deepsig_chain']}.\n"
        f"- Concordance: {nterm['n_terminal_evidence']}; final start={nterm['construct_start']}.\n"
        "- N-terminal hydropathy: hydrophobic transport-label core followed by a polar transition around residues 24-28.\n\n"
        "UniProt annotation 是数据库专家结合实验、序列和文献整理的人工审阅加工注释。DeepSig 是识别蛋白 N 端‘分泌运输标签’的深度学习网络，在本项目中只用于独立验证 signal region，不是 SignalP。\n\n"
        "## 4. C-terminal hydropathy\n\n"
        f"All primary endpoints used complete symmetric 7/9/11/15/21-aa windows. The primary stable band is {decision.get('primary_hydrophobic_boundary_band')} with peak {decision.get('primary_peak_endpoint')} and peak HydroConsensus {decision.get('primary_peak_consensus'):.6f}.\n\n"
        + _markdown_table(scale_rows) + "\n\n"
        "Kyte-Doolittle 给每个氨基酸一个偏水/偏油数值，再用滑动窗口寻找性质变化。完整对称窗口要求切点两侧都有完整数据，避免靠近右端时因缺失窗口产生人为右移。\n\n"
        "## 5. 3D structural boundary\n\n"
        "Fixed flanks are e-10..e and e+1..e+10; heavy-atom contacts use <5.0 Å. HIGH requires at least two of: density >= q75, continuous helix/strand, or mean RSA <0.10. pLDDT is confidence-only.\n\n"
        + _markdown_table(primary_structural.round(4)) + "\n\n"
        "Contact density 看切点两侧在三维空间中还有多少‘搭桥’；DSSP 判断是否剪断完整 α 螺旋/β 折叠；RSA 判断局部是否深埋核心；pLDDT 只表示 AlphaFold 对该局部结构的把握。\n\n"
        "## 6. Contact threshold sensitivity\n\n"
        + _markdown_table(sensitivity_view) + "\n\n"
        f"Endpoint 160 is acceptable in {int(round(float(decision.get('endpoint_160_acceptable_fraction')) * 7))}/7 settings, therefore sensitivity is **{decision.get('structural_sensitivity')}**. This analysis changes confidence only; q75 remains the primary veto.\n\n"
        "## 7. Final engineering decision\n\n"
        "Reviewed UniProt + DeepSig + N-terminal hydropathy → start 27; complete-window multiscale hydropathy → band 154-160; primary q75 3D veto → acceptable 154,155,156,160; maximal retention → end 160.\n\n"
        "Maximal retention 的意思是：多个相邻位置都安全时选择靠后的一个，以保留更多 IL-24 抗原序列。\n\n"
        "## 8. Why not shorter?\n\n"
        f"Structurally acceptable endpoints inside the primary band: {acceptable}. Length was used only here as the final tie-break; the selected endpoint is {decision.get('recommended_end')}.\n\n"
        "154、155、156 和 160 都在同一疏水边界带且 q75 可接受；选择最右端 160 是为了保留更多抗原表面。\n\n"
        "## 9. Why not longer?\n\n161 以后无法同时满足最大 21-aa 完整对称窗口，因此属于 edge exploratory，不与 154-160 在同一证据质量下比较。\n\n"
        "## 10. Relationship to Phase 1B\n\nPhase 1A 决定宽免疫原 construct；Phase 1B 在蛋白内部寻找 hotspot。现有 Phase 1B 热点结果未修改。\n\n"
        "## 11. Limitations\n\n结构来自 AlphaFold；DeepSig 是计算预测；constructability 不等于实验 expression yield；局部 3D risk 不等于真实 folding stability。最终仍需表达、纯化和免疫实验验证。\n",
        encoding="utf-8",
    )


def run_analysis(output_root: Path | None = None) -> dict[str, object]:
    output_root = MODULE_ROOT if output_root is None else Path(output_root)
    bands = pd.read_csv(BANDS)
    risk = pd.read_csv(RISK)
    nterm = pd.read_csv(NTERM)
    sensitivity = pd.read_csv(SENSITIVITY)
    construct_start = int(nterm.iloc[0]["construct_start"])
    decision, table = select_final_boundary(bands, risk, construct_start)
    selected_sensitivity = sensitivity.loc[sensitivity["endpoint"] == decision.get("recommended_end")].iloc[0]
    decision["n_terminal_evidence"] = str(nterm.iloc[0]["n_terminal_evidence"])
    decision["endpoint_160_acceptable_fraction"] = float(selected_sensitivity["acceptable_fraction"])
    decision["structural_sensitivity"] = str(selected_sensitivity["sensitivity_level"])
    table = table.merge(
        sensitivity[["endpoint", "acceptable_fraction", "sensitivity_level"]],
        on="endpoint",
        how="left",
        validate="one_to_one",
    )
    result_dir = output_root / "results/04_final_construct"
    result_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(result_dir / "candidate_boundary_decision.csv", index=False)
    (result_dir / "recommended_construct.json").write_text(
        json.dumps(decision, indent=2) + "\n", encoding="utf-8"
    )
    _save_figure(table, decision, output_root / "figures/figure4_final_boundary_decision.png")
    _write_report(decision, table, output_root / "reports/phase1a_construct_v2_report.md")
    return decision


def build_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(description=__doc__)


def main() -> None:
    build_parser().parse_args()
    print(json.dumps(run_analysis(), sort_keys=True))


if __name__ == "__main__":
    main()

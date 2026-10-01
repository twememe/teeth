from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd


MODULE_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = MODULE_ROOT / "scripts/final_decision.py"


def load_module():
    assert SCRIPT.is_file(), "final_decision.py is not implemented"
    spec = importlib.util.spec_from_file_location("v2_final_decision", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sample_bands() -> pd.DataFrame:
    return pd.DataFrame(
        [{"hydropathy_rank": 1, "is_primary_band": True, "band_id": "Band_A",
          "band_start": 100, "band_end": 104, "peak_endpoint": 102,
          "peak_consensus": 0.95, "mean_consensus": 0.88}]
    )


def sample_risk() -> pd.DataFrame:
    return pd.DataFrame(
        {"endpoint": [100, 101, 102, 103, 104],
         "structural_risk": ["ACCEPTABLE_STRUCTURAL_RISK", "HIGH_STRUCTURAL_RISK",
                             "ACCEPTABLE_STRUCTURAL_RISK", "HIGH_STRUCTURAL_RISK",
                             "ACCEPTABLE_STRUCTURAL_RISK"]}
    )


def test_length_is_used_only_as_final_within_band_tiebreak():
    module = load_module()
    decision, table = module.select_final_boundary(sample_bands(), sample_risk(), construct_start=27)
    assert decision["status"] == "RECOMMENDED"
    assert decision["recommended_end"] == 104
    assert decision["recommended_construct"] == "27-104"
    assert table.loc[table["is_structurally_acceptable"], "endpoint"].tolist() == [100, 102, 104]


def test_no_stable_band_returns_explicit_status_without_substitution():
    module = load_module()
    decision, table = module.select_final_boundary(pd.DataFrame(), sample_risk(), 27)
    assert decision["status"] == "NO_STABLE_C_TERMINAL_BOUNDARY"
    assert decision["recommended_end"] is None
    assert table.empty


def test_no_fallback2():
    module = load_module()
    decision, _ = module.select_final_boundary(pd.DataFrame(), sample_risk(), 27)
    assert decision == {
        "status": "NO_STABLE_C_TERMINAL_BOUNDARY",
        "start": 27,
        "recommended_end": None,
        "recommended_construct": None,
    }


def test_all_high_risk_returns_structural_conflict_without_cross_band_choice():
    module = load_module()
    risk = sample_risk()
    risk["structural_risk"] = "HIGH_STRUCTURAL_RISK"
    decision, _ = module.select_final_boundary(sample_bands(), risk, 27)
    assert decision["status"] == "STRUCTURAL_CONFLICT"
    assert decision["recommended_end"] is None


def test_formal_result_is_primary_band_154_160_and_end_160(tmp_path):
    module = load_module()
    result = module.run_analysis(output_root=tmp_path)
    assert result["start"] == 27
    assert result["primary_hydrophobic_boundary_band"] == [154, 160]
    assert result["acceptable_endpoints"] == [154, 155, 156, 160]
    assert result["recommended_end"] == 160
    assert result["recommended_construct"] == "27-160"
    assert result["construct_length_aa"] == 134
    for path in (
        "results/04_final_construct/candidate_boundary_decision.csv",
        "results/04_final_construct/recommended_construct.json",
        "figures/figure4_final_boundary_decision.png",
        "reports/phase1a_construct_v2_report.md",
    ):
        assert (tmp_path / path).is_file()
    saved = json.loads((tmp_path / "results/04_final_construct/recommended_construct.json").read_text())
    assert saved["legacy_global_selector_used"] is False
    assert saved["n_terminal_evidence"] == "STRONGLY_CONCORDANT"
    assert saved["endpoint_160_acceptable_fraction"] >= 0
    report = (tmp_path / "reports/phase1a_construct_v2_report.md").read_text(encoding="utf-8")
    assert "Phase 1A v2.1" in report
    assert "UNAVAILABLE_EXTERNAL_LICENSED_PACKAGE" in report
    assert "DeepSig" in report
    assert "Contact threshold sensitivity" in report
    assert "Phase 1B" in report


def test_new_production_scripts_have_zero_legacy_identifiers():
    forbidden = [
        "Par" + "eto",
        "Length" + "RetentionScore",
        "fallback" + "_2",
        "construct_candidates" + "_pre_evaluation",
        "pareto" + "_endpoints",
    ]
    for script in (MODULE_ROOT / "scripts").glob("*.py"):
        text = script.read_text(encoding="utf-8")
        assert all(token not in text for token in forbidden), f"legacy identifier in {script.name}"

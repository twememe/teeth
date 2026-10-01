from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


MODULE_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = MODULE_ROOT / "scripts/structural_boundary.py"


def load_module():
    assert SCRIPT.is_file(), "structural_boundary.py is not implemented"
    spec = importlib.util.spec_from_file_location("v2_structural_boundary", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_fixed_local_structure_flanks():
    module = load_module()
    for endpoint in (50, 100, 160):
        left, right = module.local_flanks(endpoint)
        assert left == list(range(endpoint - 10, endpoint + 1))
        assert right == list(range(endpoint + 1, endpoint + 11))
        assert len(left) == 11 and len(right) == 10


def test_contact_density_uses_fixed_comparable_pair_denominator():
    module = load_module()
    atoms = {position: [np.array([position * 100.0, 0.0, 0.0])] for position in range(1, 31)}
    atoms[16] = [np.array([0.0, 0.0, 0.0])]
    atoms[17] = [np.array([0.0, 0.0, 0.0])]
    atoms[21] = [np.array([1.0, 0.0, 0.0])]
    result = module.local_contact_metrics(atoms, endpoint=20)
    assert result["local_cross_contact_count"] == 1  # 16-21; 17-21 is locally excluded.
    assert result["comparable_residue_pair_count"] == module.comparable_pair_count()
    assert result["local_cross_contact_density"] == 1 / module.comparable_pair_count()
    assert module.comparable_pair_count() == module.comparable_pair_count(endpoint=150)


def test_dssp_boundary_status_is_categorical_not_scored():
    module = load_module()
    assert module.dssp_boundary_status("C", "C") == "COIL_BOUNDARY"
    assert module.dssp_boundary_status("H", "G") == "CONTINUOUS_HELIX"
    assert module.dssp_boundary_status("E", "E") == "CONTINUOUS_STRAND"
    assert module.dssp_boundary_status("H", "C") == "STRUCTURE_TRANSITION"


def test_rsa_and_plddt_are_annotations():
    module = load_module()
    assert module.rsa_boundary_status(0.099) == "DEEPLY_BURIED_BOUNDARY"
    assert module.rsa_boundary_status(0.10) == "NOT_DEEPLY_BURIED_BOUNDARY"
    assert module.structure_confidence(69.9) == "LOW_STRUCTURE_CONFIDENCE"
    assert module.structure_confidence(70.0) == "STANDARD_STRUCTURE_CONFIDENCE"


def test_high_risk_requires_two_of_three_conditions():
    module = load_module()
    frame = pd.DataFrame(
        {
            "endpoint": [10, 11, 12, 13],
            "local_cross_contact_density": [0.1, 0.2, 0.3, 0.4],
            "dssp_status": ["COIL_BOUNDARY", "CONTINUOUS_HELIX", "COIL_BOUNDARY", "CONTINUOUS_STRAND"],
            "boundary_mean_RSA": [0.5, 0.5, 0.05, 0.5],
            "boundary_median_pLDDT": [90.0, 90.0, 60.0, 60.0],
        }
    )
    result = module.assign_structural_risk(frame)
    assert result.set_index("endpoint").loc[13, "structural_risk"] == "HIGH_STRUCTURAL_RISK"
    assert result.set_index("endpoint").loc[12, "structural_risk"] == "ACCEPTABLE_STRUCTURAL_RISK"
    assert result.set_index("endpoint").loc[12, "structure_confidence"] == "LOW_STRUCTURE_CONFIDENCE"
    assert result.set_index("endpoint").loc[13, "risk_condition_count"] == 2


def test_only_hydropathy_band_endpoints_are_analyzed():
    module = load_module()
    bands = pd.DataFrame(
        {"band_start": [50, 80], "band_end": [52, 81], "hydropathy_rank": [1, 2]}
    )
    assert module.candidate_endpoints(bands) == [50, 51, 52, 80, 81]


def test_formal_outputs_exist_and_use_candidate_endpoints(tmp_path):
    module = load_module()
    result = module.run_analysis(output_root=tmp_path)
    assert result["candidate_endpoint_count"] > 0
    for path in (
        "results/03_structural_boundary/local_boundary_contacts.csv",
        "results/03_structural_boundary/dssp_boundary_status.csv",
        "results/03_structural_boundary/rsa_boundary_status.csv",
        "results/03_structural_boundary/structural_risk.csv",
        "figures/figure3_structural_risk.png",
        "figures/04_contact_threshold_sensitivity.png",
    ):
        assert (tmp_path / path).is_file()


def test_contact_quantile_sensitivity_uses_exact_requested_quantiles():
    module = load_module()
    frame = pd.DataFrame(
        {
            "endpoint": list(range(150, 161)),
            "local_cross_contact_density": [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11],
            "dssp_status": ["COIL_BOUNDARY"] * 7 + ["CONTINUOUS_HELIX"] * 4,
            "boundary_mean_RSA": [0.5] * 8 + [0.05] * 3,
            "boundary_median_pLDDT": [95.0] * 11,
        }
    )

    result = module.contact_quantile_sensitivity(frame, endpoints=range(154, 161))

    assert result["endpoint"].tolist() == list(range(154, 161))
    assert [column for column in result if column.startswith("q") and column[1:].isdigit()] == [
        "q60", "q65", "q70", "q75", "q80", "q85", "q90"
    ]
    assert result["acceptable_fraction"].between(0.0, 1.0).all()


def test_contact_quantile_sensitivity_recounts_two_of_three_and_keeps_plddt_out():
    module = load_module()
    frame = pd.DataFrame(
        {
            "endpoint": [154, 155, 156, 157, 158, 159, 160],
            "local_cross_contact_density": [0.01, 0.02, 0.03, 0.08, 0.09, 0.10, 0.04],
            "dssp_status": ["COIL_BOUNDARY", "COIL_BOUNDARY", "COIL_BOUNDARY", "CONTINUOUS_HELIX", "CONTINUOUS_HELIX", "CONTINUOUS_HELIX", "COIL_BOUNDARY"],
            "boundary_mean_RSA": [0.5, 0.5, 0.5, 0.5, 0.05, 0.05, 0.5],
            "boundary_median_pLDDT": [1.0, 99.0, 1.0, 99.0, 1.0, 99.0, 1.0],
        }
    )

    result = module.contact_quantile_sensitivity(frame).set_index("endpoint")

    assert result.loc[154, "q60"] == "ACCEPTABLE_STRUCTURAL_RISK"
    assert result.loc[158, "q90"] == "HIGH_STRUCTURAL_RISK"
    assert result.loc[154, "acceptable_fraction"] == 1.0
    assert result.loc[158, "acceptable_fraction"] == 0.0


def test_sensitivity_q75_reproduces_primary_structural_rule():
    module = load_module()
    frame = pd.DataFrame(
        {
            "endpoint": list(range(154, 161)),
            "local_cross_contact_density": [0.01, 0.02, 0.03, 0.08, 0.09, 0.10, 0.04],
            "dssp_status": ["COIL_BOUNDARY", "COIL_BOUNDARY", "COIL_BOUNDARY", "CONTINUOUS_HELIX", "CONTINUOUS_HELIX", "CONTINUOUS_HELIX", "COIL_BOUNDARY"],
            "boundary_mean_RSA": [0.5, 0.5, 0.5, 0.5, 0.05, 0.05, 0.5],
            "boundary_median_pLDDT": [95.0] * 7,
        }
    )
    primary = module.assign_structural_risk(frame).set_index("endpoint")["structural_risk"]
    sensitivity = module.contact_quantile_sensitivity(frame).set_index("endpoint")["q75"]
    assert sensitivity.to_dict() == primary.to_dict()


def test_formal_run_writes_contact_quantile_sensitivity(tmp_path):
    module = load_module()
    module.run_analysis(output_root=tmp_path)
    assert (tmp_path / "results/03_structural_boundary/contact_quantile_sensitivity.csv").is_file()
    table = pd.read_csv(tmp_path / "results/03_structural_boundary/contact_quantile_sensitivity.csv")
    assert table["endpoint"].tolist() == list(range(154, 161))

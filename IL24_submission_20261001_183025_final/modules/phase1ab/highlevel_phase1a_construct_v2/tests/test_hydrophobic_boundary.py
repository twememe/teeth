from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


MODULE_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = MODULE_ROOT / "scripts/hydrophobic_boundary.py"


def load_module():
    assert SCRIPT.is_file(), "hydrophobic_boundary.py is not implemented"
    spec = importlib.util.spec_from_file_location("v2_hydrophobic_boundary", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_complete_windows_only():
    module = load_module()
    sequence = module.read_fasta(MODULE_ROOT / "inputs/reference_181.fasta")
    table = module.build_endpoint_table(sequence, construct_start=27)
    assert table["endpoint"].tolist() == list(range(47, 161))
    assert (table["complete_scale_count"] == 5).all()
    assert table[[f"DeltaH_{w}" for w in (7, 9, 11, 15, 21)]].notna().all().all()


def test_primary_max_endpoint_is_derived_from_length_and_largest_window():
    module = load_module()
    assert module.primary_endpoint_bounds(181, 27, (7, 9, 11, 15, 21)) == (47, 160)
    assert module.primary_endpoint_bounds(100, 11, (3, 5, 9)) == (19, 91)


def test_delta_is_after_minus_before_with_equal_windows():
    module = load_module()
    values = np.array([-2.0, -2.0, 3.0, 5.0])
    assert module.delta_h(values, endpoint=2, window=2) == 6.0


def test_peak_and_band_rules_do_not_cross_valleys():
    module = load_module()
    table = pd.DataFrame(
        {
            "endpoint": list(range(10, 20)),
            "HydroConsensus": [0.4, 0.76, 0.82, 0.78, 0.5, 0.4, 0.77, 0.91, 0.8, 0.4],
            "HydroSupportCount75": [0, 3, 5, 3, 1, 0, 3, 5, 4, 0],
            "HydroSupportCount90": [0, 0, 2, 0, 0, 0, 0, 4, 1, 0],
        }
    )
    peaks, bands = module.find_peaks_and_bands(table)
    assert peaks["endpoint"].tolist() == [12, 17]
    assert bands[["band_start", "band_end"]].values.tolist() == [[11, 13], [16, 18]]


def test_no_variable_tail_statistics_in_ranking():
    module = load_module()
    sequence = module.read_fasta(MODULE_ROOT / "inputs/reference_181.fasta")
    table = module.build_endpoint_table(sequence, construct_start=27)
    forbidden = ("tail", "length", "coverage", "removed")
    assert not any(any(token in column.lower() for token in forbidden) for column in table.columns)


def test_length_not_used_before_final_tiebreak():
    module = load_module()
    bands = pd.DataFrame(
        [
            {"band_id": "A", "band_start": 50, "band_end": 55, "peak_consensus": 0.90,
             "mean_consensus": 0.82, "peak_support75": 5, "peak_support90": 3, "band_width": 6},
            {"band_id": "B", "band_start": 140, "band_end": 150, "peak_consensus": 0.89,
             "mean_consensus": 0.84, "peak_support75": 5, "peak_support90": 4, "band_width": 11},
        ]
    )
    ranked = module.rank_bands(bands)
    assert ranked.iloc[0]["band_id"] == "A"
    shifted = bands.copy()
    shifted[["band_start", "band_end"]] = shifted[["band_start", "band_end"]].iloc[::-1].to_numpy()
    assert module.rank_bands(shifted)["band_id"].tolist() == ranked["band_id"].tolist()


def test_synthetic_middle_transition_does_not_right_shift():
    module = load_module()
    sequence = "D" * 75 + "I" * 30 + "D" * 76
    table = module.build_endpoint_table(sequence, construct_start=27)
    peaks, bands = module.find_peaks_and_bands(table)
    assert not peaks.empty
    primary = module.rank_bands(bands).iloc[0]
    assert 70 <= int(primary["peak_endpoint"]) <= 80
    assert int(primary["band_end"]) < 120


def test_no_legacy_fallback_identifier_and_required_outputs(tmp_path):
    module = load_module()
    source = SCRIPT.read_text(encoding="utf-8").lower()
    assert "fallback_2" not in source
    result = module.run_analysis(output_root=tmp_path)
    assert result["status"] in {"BANDS_FOUND", "NO_STABLE_C_TERMINAL_BOUNDARY"}
    for path in (
        "results/02_hydrophobic_boundary/hydropathy_multiscale.csv",
        "results/02_hydrophobic_boundary/hydropathy_endpoint_complete_windows.csv",
        "results/02_hydrophobic_boundary/hydrophobic_peaks.csv",
        "results/02_hydrophobic_boundary/hydrophobic_boundary_bands.csv",
        "figures/figure1_hydropathy_multiscale.png",
        "figures/figure2_hydro_consensus.png",
    ):
        assert (tmp_path / path).is_file()


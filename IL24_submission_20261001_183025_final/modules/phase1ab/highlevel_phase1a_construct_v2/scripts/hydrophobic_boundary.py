"""V2 complete-window hydrophobic transition and stable-band discovery."""

from __future__ import annotations

import argparse
import csv
from functools import cmp_to_key
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


MODULE_ROOT = Path(__file__).resolve().parents[1]
FASTA = MODULE_ROOT / "inputs/reference_181.fasta"
WINDOWS = (7, 9, 11, 15, 21)
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def read_fasta(path: Path = FASTA) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    sequence = "".join(line.strip() for line in lines if line.strip() and not line.startswith(">"))
    if not sequence or set(sequence) - set(KD):
        raise ValueError("invalid protein FASTA")
    return sequence


def kd_values(sequence: str) -> np.ndarray:
    return np.array([KD[aa] for aa in sequence], dtype=float)


def primary_endpoint_bounds(
    sequence_length: int, construct_start: int, windows: Iterable[int] = WINDOWS
) -> tuple[int, int]:
    largest = max(windows)
    return construct_start + largest - 1, sequence_length - largest


def delta_h(values: np.ndarray, endpoint: int, window: int) -> float:
    before = values[endpoint - window : endpoint]
    after = values[endpoint : endpoint + window]
    if len(before) != window or len(after) != window:
        raise ValueError("incomplete symmetric window")
    return float(np.mean(after) - np.mean(before))


def inclusive_percentile(values: pd.Series) -> pd.Series:
    rounded = values.round(12)
    array = rounded.to_numpy(dtype=float)
    return pd.Series([(array <= value).sum() / len(array) for value in array], index=values.index)


def centered_profile(values: np.ndarray, window: int) -> np.ndarray:
    result = np.full(len(values), np.nan)
    half = window // 2
    for center in range(half, len(values) - half):
        result[center] = float(np.mean(values[center - half : center + half + 1]))
    return result


def build_residue_profile(sequence: str) -> pd.DataFrame:
    values = kd_values(sequence)
    data: dict[str, object] = {
        "position": np.arange(1, len(sequence) + 1),
        "AA": list(sequence),
    }
    for window in WINDOWS:
        data[f"KD_{window}"] = centered_profile(values, window)
    return pd.DataFrame(data)


def build_endpoint_table(sequence: str, construct_start: int) -> pd.DataFrame:
    values = kd_values(sequence)
    first, last = primary_endpoint_bounds(len(sequence), construct_start, WINDOWS)
    rows = []
    for endpoint in range(first, last + 1):
        row: dict[str, object] = {"endpoint": endpoint}
        for window in WINDOWS:
            row[f"DeltaH_{window}"] = delta_h(values, endpoint, window)
        row["complete_scale_count"] = len(WINDOWS)
        rows.append(row)
    table = pd.DataFrame(rows)
    percentile_columns = []
    for window in WINDOWS:
        column = f"P{window}"
        table[column] = inclusive_percentile(table[f"DeltaH_{window}"])
        percentile_columns.append(column)
    table["HydroConsensus"] = table[percentile_columns].median(axis=1)
    table["HydroSupportCount75"] = (table[percentile_columns] >= 0.75).sum(axis=1)
    table["HydroSupportCount90"] = (table[percentile_columns] >= 0.90).sum(axis=1)
    return table


def _qualifying_local_peaks(table: pd.DataFrame) -> pd.DataFrame:
    candidates = []
    rows = table.reset_index(drop=True)
    for index in range(1, len(rows) - 1):
        current = float(rows.loc[index, "HydroConsensus"])
        left = float(rows.loc[index - 1, "HydroConsensus"])
        right = float(rows.loc[index + 1, "HydroConsensus"])
        local = current >= left and current >= right and (current > left or current > right)
        strong = current >= 0.75 and int(rows.loc[index, "HydroSupportCount75"]) >= 4
        if local and strong:
            candidates.append(rows.loc[index].to_dict())
    if not candidates:
        return pd.DataFrame(columns=list(table.columns))
    peaks = pd.DataFrame(candidates)
    # Collapse adjacent equal-consensus plateau representatives without rightward preference.
    kept = []
    for _, group in peaks.groupby((peaks["endpoint"].diff().fillna(2) > 1).cumsum()):
        ordered = group.sort_values(
            ["HydroConsensus", "HydroSupportCount90", "HydroSupportCount75", "endpoint"],
            ascending=[False, False, False, True],
        )
        kept.append(ordered.iloc[0])
    return pd.DataFrame(kept).sort_values("endpoint").reset_index(drop=True)


def find_peaks_and_bands(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    peaks = _qualifying_local_peaks(table)
    if peaks.empty:
        return peaks, pd.DataFrame(
            columns=["band_id", "band_start", "band_end", "peak_endpoint", "peak_consensus",
                     "mean_consensus", "peak_support75", "peak_support90", "band_width"]
        )
    by_endpoint = table.set_index("endpoint")
    eligible = (
        (table["HydroConsensus"] >= 0.75)
        & (table["HydroSupportCount75"] >= 3)
    )
    eligible_map = dict(zip(table["endpoint"].astype(int), eligible))
    first = int(table["endpoint"].min())
    last = int(table["endpoint"].max())
    band_ranges: dict[tuple[int, int], list[int]] = {}
    for peak in peaks["endpoint"].astype(int):
        left = peak
        right = peak
        while left - 1 >= first and eligible_map.get(left - 1, False):
            left -= 1
        while right + 1 <= last and eligible_map.get(right + 1, False):
            right += 1
        band_ranges.setdefault((left, right), []).append(peak)
    rows = []
    for band_number, ((left, right), band_peaks) in enumerate(sorted(band_ranges.items()), start=1):
        members = by_endpoint.loc[left:right]
        peak_members = peaks.loc[peaks["endpoint"].isin(band_peaks)].sort_values(
            ["HydroConsensus", "HydroSupportCount90", "endpoint"],
            ascending=[False, False, True],
        )
        peak = peak_members.iloc[0]
        rows.append(
            {
                "band_id": f"Band_{band_number}",
                "band_start": left,
                "band_end": right,
                "peak_endpoint": int(peak["endpoint"]),
                "peak_consensus": float(peak["HydroConsensus"]),
                "mean_consensus": float(members["HydroConsensus"].mean()),
                "peak_support75": int(peak["HydroSupportCount75"]),
                "peak_support90": int(peak["HydroSupportCount90"]),
                "band_width": right - left + 1,
            }
        )
    return peaks, pd.DataFrame(rows)


def _compare_bands(left: dict, right: dict) -> int:
    peak_difference = float(left["peak_consensus"]) - float(right["peak_consensus"])
    if abs(peak_difference) > 0.01:
        return -1 if peak_difference > 0 else 1
    keys = (
        ("mean_consensus", True),
        ("peak_support90", True),
        ("peak_support75", True),
        ("band_width", False),
    )
    for key, larger_is_better in keys:
        a, b = float(left[key]), float(right[key])
        if a != b:
            better = a > b if larger_is_better else a < b
            return -1 if better else 1
    return -1 if str(left["band_id"]) < str(right["band_id"]) else 1


def rank_bands(bands: pd.DataFrame) -> pd.DataFrame:
    if bands.empty:
        result = bands.copy()
        result["hydropathy_rank"] = pd.Series(dtype=int)
        result["is_primary_band"] = pd.Series(dtype=bool)
        return result
    records = sorted(bands.to_dict("records"), key=cmp_to_key(_compare_bands))
    result = pd.DataFrame(records)
    result.insert(0, "hydropathy_rank", range(1, len(result) + 1))
    result.insert(1, "is_primary_band", result["hydropathy_rank"] == 1)
    return result


def build_edge_exploratory(sequence: str) -> pd.DataFrame:
    values = kd_values(sequence)
    last_primary = len(sequence) - max(WINDOWS)
    rows = []
    for endpoint in range(last_primary + 1, len(sequence)):
        row: dict[str, object] = {"endpoint": endpoint, "stage": "edge_exploratory"}
        complete = 0
        for window in WINDOWS:
            if endpoint + window <= len(sequence):
                row[f"DeltaH_{window}"] = delta_h(values, endpoint, window)
                complete += 1
            else:
                row[f"DeltaH_{window}"] = np.nan
        row["complete_scale_count"] = complete
        rows.append(row)
    return pd.DataFrame(rows)


def _save_figures(
    profile: pd.DataFrame,
    endpoints: pd.DataFrame,
    peaks: pd.DataFrame,
    bands: pd.DataFrame,
    figure_dir: Path,
) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)
    for window in WINDOWS:
        ax.plot(profile["position"], profile[f"KD_{window}"], label=f"KD {window}")
    ax.axvspan(1, 26, color="#f4a261", alpha=0.18, label="signal peptide")
    ax.axvspan(int(endpoints.endpoint.min()), int(endpoints.endpoint.max()), color="#90be6d", alpha=0.10, label="primary valid endpoints")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set(xlabel="Project residue", ylabel="Mean Kyte-Doolittle hydropathy", title="Complete-window multiscale hydropathy")
    ax.legend(ncol=4, frameon=False)
    fig.savefig(figure_dir / "figure1_hydropathy_multiscale.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 5.8), constrained_layout=True)
    ax.plot(endpoints["endpoint"], endpoints["HydroConsensus"], color="#277da1", linewidth=2)
    ax.axhline(0.75, color="#666666", linestyle="--", label="band threshold 0.75")
    if not peaks.empty:
        ax.scatter(peaks["endpoint"], peaks["HydroConsensus"], marker="^", s=70, color="#d00000", label="qualifying local peaks")
    for row in bands.itertuples():
        color = "#f9c74f" if bool(row.is_primary_band) else "#90be6d"
        ax.axvspan(row.band_start - 0.5, row.band_end + 0.5, color=color, alpha=0.25)
    ax.set(xlabel="Complete-window endpoint", ylabel="HydroConsensus", title="Stable hydrophobic boundary bands")
    ax.legend(frameon=False)
    fig.savefig(figure_dir / "figure2_hydro_consensus.png", dpi=220)
    plt.close(fig)


def run_analysis(output_root: Path | None = None) -> dict[str, object]:
    output_root = MODULE_ROOT if output_root is None else Path(output_root)
    sequence = read_fasta()
    profile = build_residue_profile(sequence)
    endpoints = build_endpoint_table(sequence, construct_start=27)
    peaks, bands = find_peaks_and_bands(endpoints)
    ranked_bands = rank_bands(bands)
    edge = build_edge_exploratory(sequence)
    result_dir = output_root / "results/02_hydrophobic_boundary"
    result_dir.mkdir(parents=True, exist_ok=True)
    profile.to_csv(result_dir / "hydropathy_multiscale.csv", index=False, na_rep="NA")
    endpoints.to_csv(result_dir / "hydropathy_endpoint_complete_windows.csv", index=False)
    edge.to_csv(result_dir / "edge_exploratory.csv", index=False, na_rep="NA")
    peaks.to_csv(result_dir / "hydrophobic_peaks.csv", index=False)
    ranked_bands.to_csv(result_dir / "hydrophobic_boundary_bands.csv", index=False)
    _save_figures(profile, endpoints, peaks, ranked_bands, output_root / "figures")
    status = "BANDS_FOUND" if not ranked_bands.empty else "NO_STABLE_C_TERMINAL_BOUNDARY"
    primary = ranked_bands.iloc[0].to_dict() if not ranked_bands.empty else None
    return {
        "status": status,
        "primary_min_endpoint": int(endpoints.endpoint.min()),
        "primary_max_endpoint": int(endpoints.endpoint.max()),
        "peak_count": len(peaks),
        "band_count": len(ranked_bands),
        "primary_band": primary,
    }


def build_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(description=__doc__)


def main() -> None:
    build_parser().parse_args()
    result = run_analysis()
    for key, value in result.items():
        print(f"{key}={value}")


if __name__ == "__main__":
    main()

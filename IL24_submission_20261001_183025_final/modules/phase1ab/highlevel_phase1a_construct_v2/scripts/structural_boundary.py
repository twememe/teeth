"""V2 fixed-flank structural boundary risk analysis."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


MODULE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = MODULE_ROOT.parent
PDB = MODULE_ROOT / "inputs/reference_structure_181.pdb"
BANDS = MODULE_ROOT / "results/02_hydrophobic_boundary/hydrophobic_boundary_bands.csv"
MKDSSP = __import__("os").environ.get("MKDSSP", "mkdssp")
WSL_PYTHON = __import__("sys").executable
SENSITIVITY_QUANTILES = (60, 65, 70, 75, 80, 85, 90)
PRIMARY_BAND_ENDPOINTS = tuple(range(154, 161))


def local_flanks(endpoint: int) -> tuple[list[int], list[int]]:
    return list(range(endpoint - 10, endpoint + 1)), list(range(endpoint + 1, endpoint + 11))


def comparable_pair_count(endpoint: int = 100, exclude_leq: int = 4) -> int:
    left, right = local_flanks(endpoint)
    return sum(1 for i in left for j in right if abs(i - j) > exclude_leq)


def local_contact_metrics(
    atoms_by_position: dict[int, list[np.ndarray]],
    endpoint: int,
    cutoff: float = 5.0,
    exclude_leq: int = 4,
) -> dict[str, object]:
    left, right = local_flanks(endpoint)
    comparable = [(i, j) for i in left for j in right if abs(i - j) > exclude_leq]
    contacts = []
    minimum_distances = []
    for i, j in comparable:
        distances = [
            float(np.linalg.norm(atom_i - atom_j))
            for atom_i in atoms_by_position[i]
            for atom_j in atoms_by_position[j]
        ]
        minimum = min(distances)
        if minimum < cutoff:
            contacts.append((i, j))
            minimum_distances.append(minimum)
    count = len(contacts)
    denominator = len(comparable)
    return {
        "endpoint": endpoint,
        "left_flank_start": left[0],
        "left_flank_end": left[-1],
        "right_flank_start": right[0],
        "right_flank_end": right[-1],
        "left_flank_residue_count": len(left),
        "right_flank_residue_count": len(right),
        "comparable_residue_pair_count": denominator,
        "local_cross_contact_count": count,
        "local_cross_contact_density": count / denominator,
        "minimum_contact_distance_A": min(minimum_distances) if minimum_distances else np.nan,
    }


def parse_pdb(path: Path = PDB) -> tuple[dict[int, list[np.ndarray]], dict[int, float]]:
    atoms: dict[int, list[np.ndarray]] = {}
    plddt: dict[int, float] = {}
    seen_residues: list[int] = []
    for line in path.read_text(encoding="ascii").splitlines():
        if not line.startswith("ATOM  ") or line[21:22] != "A":
            continue
        position = int(line[22:26])
        atom_name = line[12:16].strip()
        element = line[76:78].strip() or atom_name[0]
        if position not in atoms:
            atoms[position] = []
            seen_residues.append(position)
        if element.upper() != "H":
            atoms[position].append(
                np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            )
        if atom_name == "CA":
            plddt[position] = float(line[60:66])
    if seen_residues != list(range(1, 182)):
        raise ValueError("PDB Chain A is not continuous 1-181")
    return atoms, plddt


def _run_wsl(arguments: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        arguments if __import__("os").name != "nt" else ["wsl.exe", "-d", "Ubuntu", "--cd", str(WORKSPACE), "--", *arguments],
        cwd=WORKSPACE,
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )


def parse_dssp(text: str) -> dict[int, str]:
    labels: dict[int, str] = {}
    started = False
    for line in text.splitlines():
        if "#  RESIDUE AA STRUCTURE" in line:
            started = True
            continue
        if not started or len(line) < 17:
            continue
        try:
            position = int(line[5:10])
        except ValueError:
            continue
        if line[11:12] != "A":
            continue
        labels[position] = line[16:17].strip() or "C"
    return labels


def run_dssp() -> tuple[dict[int, str], str]:
    relative_pdb = PDB.relative_to(WORKSPACE).as_posix()
    completed = _run_wsl([MKDSSP, "--output-format", "dssp", relative_pdb])
    labels = parse_dssp(completed.stdout)
    if set(labels) != set(range(1, 182)):
        raise RuntimeError(f"DSSP mapped {len(labels)} residues, expected 181")
    return labels, "mkdssp 4.6.1"


def run_freesasa() -> tuple[dict[int, float], str]:
    relative_pdb = PDB.relative_to(WORKSPACE).as_posix()
    code = (
        "import json,sys,shutil,tempfile,freesasa;"
        "src=sys.argv[1];tmp=tempfile.NamedTemporaryFile(suffix='.pdb',delete=False).name;"
        "shutil.copyfile(src,tmp);s=freesasa.Structure(tmp);r=freesasa.calc(s);"
        "a=r.residueAreas()['A'];"
        "print(json.dumps({int(k):v.relativeTotal for k,v in a.items()}))"
    )
    completed = _run_wsl([WSL_PYTHON, "-c", code, relative_pdb])
    rsa = {int(key): float(value) for key, value in json.loads(completed.stdout).items()}
    if set(rsa) != set(range(1, 182)):
        raise RuntimeError(f"FreeSASA mapped {len(rsa)} residues, expected 181")
    return rsa, "FreeSASA 2.2.1"


def _dssp_class(label: str) -> str:
    if label in {"H", "G", "I"}:
        return "H"
    if label in {"E", "B"}:
        return "E"
    return "C"


def dssp_boundary_status(left_label: str, right_label: str) -> str:
    left = _dssp_class(left_label)
    right = _dssp_class(right_label)
    if left == right == "C":
        return "COIL_BOUNDARY"
    if left == right == "H":
        return "CONTINUOUS_HELIX"
    if left == right == "E":
        return "CONTINUOUS_STRAND"
    return "STRUCTURE_TRANSITION"


def rsa_boundary_status(boundary_mean: float) -> str:
    return "DEEPLY_BURIED_BOUNDARY" if boundary_mean < 0.10 else "NOT_DEEPLY_BURIED_BOUNDARY"


def structure_confidence(boundary_median_plddt: float) -> str:
    return "LOW_STRUCTURE_CONFIDENCE" if boundary_median_plddt < 70 else "STANDARD_STRUCTURE_CONFIDENCE"


def candidate_endpoints(bands: pd.DataFrame) -> list[int]:
    endpoints = {
        endpoint
        for row in bands.itertuples()
        for endpoint in range(int(row.band_start), int(row.band_end) + 1)
    }
    return sorted(endpoints)


def assign_structural_risk(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    threshold = float(result["local_cross_contact_density"].quantile(0.75))
    result["contact_density_q75"] = threshold
    result["condition_high_contact_density"] = result["local_cross_contact_density"] >= threshold
    result["condition_continuous_structure"] = result["dssp_status"].isin(
        ["CONTINUOUS_HELIX", "CONTINUOUS_STRAND"]
    )
    result["condition_deeply_buried"] = result["boundary_mean_RSA"] < 0.10
    condition_columns = [
        "condition_high_contact_density",
        "condition_continuous_structure",
        "condition_deeply_buried",
    ]
    result["risk_condition_count"] = result[condition_columns].sum(axis=1)
    result["structural_risk"] = np.where(
        result["risk_condition_count"] >= 2,
        "HIGH_STRUCTURAL_RISK",
        "ACCEPTABLE_STRUCTURAL_RISK",
    )
    result["structure_confidence"] = result["boundary_median_pLDDT"].map(structure_confidence)
    return result


def contact_quantile_sensitivity(
    frame: pd.DataFrame,
    endpoints: object = PRIMARY_BAND_ENDPOINTS,
) -> pd.DataFrame:
    endpoint_set = {int(endpoint) for endpoint in endpoints}
    selected = frame.loc[frame["endpoint"].astype(int).isin(endpoint_set)].copy()
    if sorted(selected["endpoint"].astype(int).tolist()) != sorted(endpoint_set):
        raise ValueError("sensitivity requires every requested endpoint exactly once")
    selected = selected.sort_values("endpoint").reset_index(drop=True)
    selected = selected.rename(
        columns={
            "local_cross_contact_density": "density",
            "dssp_status": "DSSP",
            "boundary_mean_RSA": "RSA",
        }
    )
    continuous = selected["DSSP"].isin(["CONTINUOUS_HELIX", "CONTINUOUS_STRAND"])
    buried = selected["RSA"] < 0.10
    for quantile in SENSITIVITY_QUANTILES:
        threshold = float(frame["local_cross_contact_density"].quantile(quantile / 100))
        selected[f"q{quantile}_threshold"] = threshold
        risk_count = (selected["density"] >= threshold).astype(int) + continuous.astype(int) + buried.astype(int)
        selected[f"q{quantile}"] = np.where(
            risk_count >= 2,
            "HIGH_STRUCTURAL_RISK",
            "ACCEPTABLE_STRUCTURAL_RISK",
        )
    risk_columns = [f"q{quantile}" for quantile in SENSITIVITY_QUANTILES]
    selected["acceptable_count"] = selected[risk_columns].eq("ACCEPTABLE_STRUCTURAL_RISK").sum(axis=1)
    selected["acceptable_fraction"] = selected["acceptable_count"] / len(risk_columns)
    selected["sensitivity_level"] = np.select(
        [selected["acceptable_count"] >= 5, selected["acceptable_count"] >= 3],
        ["ROBUST", "MODERATE"],
        default="FRAGILE",
    )
    ordered = ["endpoint", "density", "DSSP", "RSA", "boundary_median_pLDDT"]
    ordered += [item for quantile in SENSITIVITY_QUANTILES for item in (f"q{quantile}_threshold", f"q{quantile}")]
    ordered += ["acceptable_count", "acceptable_fraction", "sensitivity_level"]
    return selected[ordered]


def _save_figure(risk: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    x = risk["endpoint"]
    fig, axes = plt.subplots(4, 1, figsize=(12, 11), sharex=True, constrained_layout=True)
    axes[0].plot(x, risk["local_cross_contact_density"], marker="o")
    axes[0].axhline(risk["contact_density_q75"].iloc[0], linestyle="--", color="#d00000")
    axes[0].set_ylabel("Contact density")
    dssp_codes = {"COIL_BOUNDARY": 0, "STRUCTURE_TRANSITION": 1, "CONTINUOUS_HELIX": 2, "CONTINUOUS_STRAND": 3}
    axes[1].scatter(x, risk["dssp_status"].map(dssp_codes), c=risk["risk_condition_count"], cmap="Reds")
    axes[1].set_yticks(list(dssp_codes.values()), list(dssp_codes.keys()))
    axes[1].set_ylabel("DSSP status")
    axes[2].plot(x, risk["boundary_mean_RSA"], marker="o", color="#6a4c93")
    axes[2].axhline(0.10, linestyle="--", color="#d00000")
    axes[2].set_ylabel("Boundary mean RSA")
    axes[3].plot(x, risk["boundary_median_pLDDT"], marker="o", color="#2a9d8f")
    axes[3].axhline(70, linestyle="--", color="#666666")
    axes[3].set(ylabel="Boundary pLDDT", xlabel="Hydropathy-band endpoint")
    fig.suptitle("Fixed-flank structural risk for hydropathy candidate endpoints")
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _save_sensitivity_figure(sensitivity: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    risk = sensitivity[[f"q{q}" for q in SENSITIVITY_QUANTILES]].eq("HIGH_STRUCTURAL_RISK").astype(int)
    fig, ax = plt.subplots(figsize=(10, 5.5), constrained_layout=True)
    image = ax.imshow(risk.to_numpy(), cmap="RdYlGn_r", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(SENSITIVITY_QUANTILES)), [f"q{q}" for q in SENSITIVITY_QUANTILES])
    ax.set_yticks(range(len(sensitivity)), sensitivity["endpoint"].astype(str))
    ax.set(xlabel="Contact-density quantile", ylabel="Endpoint", title="Local contact-threshold sensitivity (green=acceptable, red=high risk)")
    for i in range(len(sensitivity)):
        for j in range(len(SENSITIVITY_QUANTILES)):
            ax.text(j, i, "HIGH" if risk.iloc[i, j] else "OK", ha="center", va="center", fontsize=8)
    fig.colorbar(image, ax=ax, ticks=[0, 1], label="Risk class")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def run_analysis(output_root: Path | None = None) -> dict[str, object]:
    output_root = MODULE_ROOT if output_root is None else Path(output_root)
    bands = pd.read_csv(BANDS)
    endpoints = candidate_endpoints(bands)
    if not endpoints:
        raise RuntimeError("no hydropathy candidate endpoints")
    atoms, plddt = parse_pdb()
    dssp, dssp_version = run_dssp()
    rsa, freesasa_version = run_freesasa()
    rows = []
    for endpoint in endpoints:
        contact = local_contact_metrics(atoms, endpoint)
        window = list(range(endpoint - 2, endpoint + 3))
        row = {
            **contact,
            "dssp_e_raw": dssp[endpoint],
            "dssp_e_plus_1_raw": dssp[endpoint + 1],
            "dssp_status": dssp_boundary_status(dssp[endpoint], dssp[endpoint + 1]),
            "RSA_e": rsa[endpoint],
            "RSA_e_plus_1": rsa[endpoint + 1],
            "boundary_mean_RSA": float(np.mean([rsa[position] for position in window])),
            "boundary_median_pLDDT": float(np.median([plddt[position] for position in window])),
        }
        row["rsa_status"] = rsa_boundary_status(float(row["boundary_mean_RSA"]))
        rows.append(row)
    risk = assign_structural_risk(pd.DataFrame(rows))
    result_dir = output_root / "results/03_structural_boundary"
    result_dir.mkdir(parents=True, exist_ok=True)
    contact_columns = [
        "endpoint", "left_flank_start", "left_flank_end", "right_flank_start", "right_flank_end",
        "left_flank_residue_count", "right_flank_residue_count", "comparable_residue_pair_count",
        "local_cross_contact_count", "local_cross_contact_density", "minimum_contact_distance_A",
    ]
    risk[contact_columns].to_csv(result_dir / "local_boundary_contacts.csv", index=False)
    risk[["endpoint", "dssp_e_raw", "dssp_e_plus_1_raw", "dssp_status"]].to_csv(
        result_dir / "dssp_boundary_status.csv", index=False
    )
    risk[["endpoint", "RSA_e", "RSA_e_plus_1", "boundary_mean_RSA", "rsa_status"]].to_csv(
        result_dir / "rsa_boundary_status.csv", index=False
    )
    risk.to_csv(result_dir / "structural_risk.csv", index=False)
    sensitivity = contact_quantile_sensitivity(risk)
    sensitivity.to_csv(result_dir / "contact_quantile_sensitivity.csv", index=False)
    _save_sensitivity_figure(sensitivity, output_root / "figures/04_contact_threshold_sensitivity.png")
    _save_figure(risk, output_root / "figures/figure3_structural_risk.png")
    (output_root / "logs").mkdir(parents=True, exist_ok=True)
    (output_root / "logs/structural_dependencies.json").write_text(
        json.dumps({"DSSP": dssp_version, "FreeSASA": freesasa_version}, indent=2) + "\n",
        encoding="utf-8",
    )
    return {
        "candidate_endpoint_count": len(endpoints),
        "high_risk_count": int((risk["structural_risk"] == "HIGH_STRUCTURAL_RISK").sum()),
        "acceptable_count": int((risk["structural_risk"] == "ACCEPTABLE_STRUCTURAL_RISK").sum()),
        "density_q75": float(risk["contact_density_q75"].iloc[0]),
    }


def build_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(description=__doc__)


def main() -> None:
    build_parser().parse_args()
    print(json.dumps(run_analysis(), sort_keys=True))


if __name__ == "__main__":
    main()

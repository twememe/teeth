#!/usr/bin/env python3
"""Incremental, read-only baseline audit and Full181 ensemble export.

No prediction, hashing, fitting, or graphics. Missing runs are never empty epitopes.
All outputs are confined to results/full181_supplement (or --output).

Server integration: original condition-specific tables precede release fallbacks;
optional raw CPU audit rows, downstream status, and frozen-prior scores stay keyed
by construct/condition/seed. Saved baseline inputs are never rewritten.
"""
import argparse
import csv
import itertools
import json
import math
import re
import shutil
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from full181_contract import CONSTRUCTS, FULL181_SEQUENCE, build_structure_mapping

NA = "NA"
CONDITIONS = ("empty", "formal_msa")
ROUTES = {"native_blind": "Native155", "immunogen_blind": "Immunogen134"}
METRICS = ("confidence_score", "ptm", "iptm", "complex_plddt",
           "epitope_residues_4p5", "paratope_residues_4p5",
           "cdr_paratope_residues_4p5", "hotspot_hits_4p5", "hotspot_coverage", "hotspot_fraction_in_interface",
           "prodigy_dG_kcal_mol", "prodigy_Kd_M_25C",
           "chain_breaks", "peptide_cn_outliers", "interchain_clashes_lt1p8",
           "nonlocal_intrachain_clashes_lt1p8")


def read(path):
    if not path.is_file():
        return []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({k: NA if row.get(k) is None else row.get(k, NA)
                         for k in fields} for row in rows)


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def residues(value):
    if value is None or str(value).strip().upper() in ("NA", "NAN", "UNKNOWN"):
        return None
    return {int(x) for x in str(value).split(";") if x.strip()}


def joined(values):
    return ";".join(map(str, sorted(values))) if values is not None else NA


def key(row):
    return row["construct"], row["condition"], int(row["seed"])


def base(construct, condition, seed, source_path):
    return dict(construct=construct, condition=condition, seed=seed,
                source_path=str(source_path))


def first_existing(*paths):
    return next((path for path in paths if path.is_file()), paths[0])


def baseline_paths(root, condition):
    original = "formal" if condition == "empty" else "formal_msa"
    paths = {"summary": root / f"results/13_interface_analysis/{original}/complex_model_summary.csv",
             "contacts": root / f"results/13_interface_analysis/{original}/interface_contacts.csv",
             "affinity": root / f"results/15_affinity/prodigy_{original}.csv"}
    if condition == "formal_msa":
        for field, flat in (("summary", "complex_model_summary.csv"),
                            ("contacts", "interface_contacts.csv"), ("affinity", "prodigy_formal_msa.csv")):
            paths[field] = first_existing(paths[field], root / "results" / flat)
    return paths


def keyed_rows(rows, label):
    result = {key(row): row for row in rows}
    if len(result) != len(rows):
        raise ValueError(f"Duplicate construct/condition/seed in {label}")
    return result


def normalized_legacy(rows, condition, source_path):
    result = []
    for saved in rows:
        construct = ROUTES[saved["route"]]
        seed = int(saved["seed"])
        if seed not in range(1, 16):
            raise ValueError(f"Unexpected legacy seed {seed} in {source_path}")
        result.append(dict(saved, construct=construct, condition=condition, seed=seed,
                           source_path=str(source_path)))
    return result


def add_prodigy(row, record, path, historical=False):
    status = record.get("status") or ("saved_per_model_table" if record and historical else "missing_observation")
    numeric_available = status in ("complete", "completed", "success", "saved_per_model_table")
    for field in ("prodigy_dG_kcal_mol", "prodigy_Kd_M_25C", "intermolecular_contacts"):
        value = number(record.get(field))
        row[field] = value if value is not None and (numeric_available or field == "intermolecular_contacts" and status == "no_interface") else NA
    row.update(prodigy_status=status,
               prodigy_failure_reason=record.get("failure_reason") or ("" if numeric_available else "PRODIGY observation unavailable"),
               prodigy_temperature_C=number(record.get("temperature_C")),
               prodigy_configured_temperature_C=number(record.get("configured_temperature_C")),
               prodigy_source_path=str(path),
               prodigy_structure_source_path=record.get("pdb_path", NA) if historical else record.get("source_path", NA),
               prodigy_evidence_level="saved_per_model_table" if record and historical else "current_original_model_score" if numeric_available else "current_prodigy_attempt" if record else "missing_observation",
               prodigy_temperature_status="reported_actual_temperature" if number(record.get("temperature_C")) is not None else "unknown_not_inferred_from_25C_column_name")


def add_geometry(row, qc, evidence):
    for field in ("chain_breaks", "peptide_cn_outliers", "interchain_clashes_lt1p8", "nonlocal_intrachain_clashes_lt1p8"):
        row[field] = qc.get(field, qc.get(field.replace("_lt1p8", "_within_1p8"), NA))
    row.update(geometry_warning=qc.get("geometry_warning", NA),
               geometry_source_path=qc.get("source_path", NA),
               geometry_evidence_level=evidence, geometry_status=qc.get("status", NA))


def stats(values):
    values = [number(v) for v in values]
    values = [v for v in values if v is not None]
    return dict(n=len(values), mean=statistics.mean(values) if values else NA,
                sd=statistics.stdev(values) if len(values) > 1 else NA,
                median=statistics.median(values) if values else NA,
                min=min(values) if values else NA, max=max(values) if values else NA)


def metadata(path):
    return dict(source_path=str(path.resolve()), exists=path.is_file(),
                size_bytes=path.stat().st_size if path.is_file() else NA,
                mtime_ns=path.stat().st_mtime_ns if path.is_file() else NA)


def pair_score(left, right):
    if left is None or right is None:
        return dict(jaccard_legacy=NA, jaccard_nonempty_union=NA,
                    interface_state="missing_observation", both_empty=NA)
    union = left | right
    value = len(left & right) / len(union) if union else 1.0
    state = "both_empty" if not union else "one_empty" if not left or not right else "both_nonempty"
    return dict(jaccard_legacy=value, jaccard_nonempty_union=value if union else NA,
                interface_state=state, both_empty=int(not union))


def load_baseline(root, audit):
    paths = {condition: baseline_paths(root, condition) for condition in CONDITIONS}
    summaries, affinity_rows = [], []
    for condition, sources in paths.items():
        current = normalized_legacy(read(sources["summary"]), condition, sources["summary"])
        if sources["summary"] == root / "results/complex_model_summary.csv":
            ranking_path = first_existing(root / "results/14_quality/model_ranking_formal_msa.csv",
                                          root / "results/model_ranking_formal_msa.csv")
            ranking = {(r["route"], int(r["seed"])): r for r in read(ranking_path)}
            # Only the release's flattened fallback needs this condition proof.
            for row in current:
                saved = ranking[(row["route"], int(row["seed"]))]
                for metric in ("confidence_score", "iptm", "complex_plddt"):
                    if abs(float(row[metric])-float(saved[metric])) > 1e-12:
                        raise ValueError("Cannot establish flattened summary condition")
        summaries.extend(current)
        affinity_rows.extend(normalized_legacy(read(sources["affinity"]), condition, sources["affinity"]))
    source = keyed_rows(summaries, "baseline summaries")
    affinity = keyed_rows(affinity_rows, "baseline PRODIGY")
    server_audit_path = audit / "server_raw_model_audit.csv"
    server_audit = keyed_rows(read(server_audit_path), "server raw-model audit")
    expected = {(c, d, s) for c,d,s in itertools.product(ROUTES.values(), CONDITIONS, range(1,16))}
    if set(server_audit) - expected:
        raise ValueError("Unexpected model key in server raw-model audit")
    pdb_rows, pdb_by_key = [], defaultdict(list)
    for model_key, checked in server_audit.items():
        path = Path(checked.get("pdb_path") or checked["source_path"])
        if not path.is_file():
            continue
        row = base(*model_key, path)
        row.update(metadata(path), evidence_level="current_original_PDB_from_server_audit",
                   validation_status=checked.get("status", NA), audit_source_path=str(server_audit_path),
                   confidence_json_available=Path(checked.get("confidence_path", NA)).is_file(),
                   pae_available=Path(checked.get("pae_path", NA)).is_file())
        pdb_rows.append(row)
        pdb_by_key[model_key].append(str(path))
    candidate_sources = {Path(r["source_pdb"]).name: r for r in read(root / "results.csv")}
    for path in sorted((root / "results").glob("*.pdb")):
        match = re.match(r"(native_blind|immunogen_blind).*seed(\d+)\.pdb$", path.name)
        if not match:
            continue
        construct, seed = ROUTES[match[1]], int(match[2])
        if "_msa_" not in path.name:
            declared = candidate_sources.get(path.name, {})
            if declared.get("msa_condition") != "formal_msa" or int(declared.get("seed", -1)) != seed:
                # Original projects may have other representative PDBs; the
                # authoritative server audit, rather than filename guesses, handles them.
                continue
        mapping = build_structure_mapping(path, construct, root)
        row = base(construct, "formal_msa", seed, path)
        row.update(metadata(path), evidence_level="current_unrelaxed_pdb_snapshot",
                   chain_lengths=";".join(f"{c}:{len(mapping[c]['residue_map'])}" for c in ("A", "H", "L")),
                   exact_sequence_validation="passed", confidence_json_available=False,
                   pae_available=False)
        pdb_rows.append(row)
        pdb_by_key[(construct, "formal_msa", seed)].append(str(path))
    unique_pdbs = {}
    for row in pdb_rows:
        existing = unique_pdbs.setdefault(row["source_path"], row)
        if key(existing) != key(row):
            raise ValueError("One original PDB path was assigned to different model keys")
    pdb_rows = list(unique_pdbs.values())
    write(audit / "available_structure_inventory.csv", pdb_rows,
          list(dict.fromkeys(k for row in pdb_rows for k in row)) or
          ["construct", "condition", "seed", "source_path", "exists"])
    models = []
    for construct, condition, seed in itertools.product(ROUTES.values(), CONDITIONS, range(1, 16)):
        row = base(construct, condition, seed, paths[condition]["summary"])
        old = source.get((construct, condition, seed))
        if old:
            row.update(old, construct=construct, condition=condition, seed=seed,
                       observation_available=True, evidence_level="saved_per_model_table",
                       status="historical_saved_model_table", missing_reason="")
        else:
            row.update(observation_available=False, evidence_level="historical_summary_only",
                       status="per_model_files_missing", missing_reason=f"No saved per-model summary for {construct}/{condition}/seed{seed} at {paths[condition]['summary']}.",
                       epitope_project_residues=NA)
        row["available_pdb_paths"] = ";".join(sorted(set(pdb_by_key[key(row)]))) or NA
        row["current_raw_prediction_bundle_verified"] = False
        add_prodigy(row, affinity.get(key(row), {}), paths[condition]["affinity"], historical=True)
        checked = server_audit.get(key(row))
        if checked:
            row.update(raw_audit_status=checked.get("status", NA), raw_audit_source_path=str(server_audit_path))
            for field in ("pdb_path", "confidence_path", "pae_path"):
                row[field] = checked.get(field, NA)
            row["current_raw_prediction_bundle_verified"] = (
                checked.get("status") in ("complete", "completed", "passed", "accepted")
                and all(Path(checked.get(field, NA)).is_file() for field in ("pdb_path", "confidence_path", "pae_path")))
            add_geometry(row, checked, "recomputed_existing_original_PDB_server_audit")
        models.append(row)
    write(audit / "baseline_per_model_inventory.csv", models)
    return models, summaries, pdb_rows


def load_full(out):
    path = out / "full181_model_summary.csv"
    rows = read(path)
    summaries = {key(r): r for r in rows}
    if len(rows) != len(summaries):
        raise ValueError("Duplicate Full181 summary key")
    expected = {("Full181", c, s) for c,s in itertools.product(CONDITIONS, range(1,16))}
    if set(summaries) - expected:
        raise ValueError("Full181 summary contains an unauthorized construct, condition or seed")
    manifest_rows = read(out / "task_manifest.csv")
    manifest = {key(r): r for r in manifest_rows}
    if len(manifest) != len(manifest_rows):
        raise ValueError("Duplicate task manifest key")
    if manifest and set(manifest) != expected:
        raise ValueError("Manifest must contain exactly the 30 authorized Full181 keys")
    result = []
    for condition, seed in itertools.product(CONDITIONS, range(1, 16)):
        task = manifest.get(("Full181", condition, seed), {})
        saved = summaries.get(("Full181", condition, seed), {})
        row = base("Full181", condition, seed, saved.get("source_path", task.get("source_path", path)))
        row.update(saved)
        row["construct"], row["condition"], row["seed"] = "Full181", condition, seed
        state = saved.get("analysis_status", saved.get("status", "missing"))
        available = state in ("completed", "complete", "success", "analyzed")
        # An explicit empty string is valid only in a completed analysis row.
        available = available and "epitope_project_residues" in saved and residues(saved["epitope_project_residues"]) is not None
        if task:
            # The analyzer accepts this state only after raw-artifact and saved
            # provenance validation; a failed downstream score must not erase it.
            available = available and task.get("status") in ("complete", "completed", "success", "downstream_blocked")
        row.update(observation_available=available,
                   evidence_level="current_original_model_analysis" if available else "missing_observation",
                   status=state if available else task.get("status", state),
                   missing_reason="" if available else saved.get("failure_reason") or task.get("failure_reason") or "No accepted Full181 model summary.",
                   task_status=task.get("status", NA), task_failure_reason=task.get("failure_reason", ""),
                   model_summary_source_path=str(path))
        if not available:
            row["epitope_project_residues"] = NA
            for metric in METRICS:
                row[metric] = NA
        result.append(row)
    affinity = keyed_rows(read(out / "full181_prodigy.csv"), "Full181 PRODIGY")
    geometry = keyed_rows(read(out / "full181_geometry_qc.csv"), "Full181 geometry")
    if (set(affinity) | set(geometry)) - expected:
        raise ValueError("Unauthorized model key in Full181 downstream tables")
    for row in result:
        add_prodigy(row, affinity.get(key(row), {}), out / "full181_prodigy.csv")
        if not row["observation_available"]:
            for field in ("prodigy_dG_kcal_mol", "prodigy_Kd_M_25C", "intermolecular_contacts"):
                row[field] = NA
        if row["observation_available"]:
            qc = geometry.get(key(row), {})
            add_geometry(row, qc, "current_original_model_analysis" if qc else "missing_observation")
    return result


def add_full_prior(models, prior, prior_path):
    for row in models:
        if row["construct"] != "Full181":
            continue
        row.update(prior_source_path=str(prior_path), prior_residue_count=len(prior) if prior else NA,
                   prior_status="available" if prior else "blocked_missing_frozen_prior")
        if prior and row["observation_available"]:
            epitope = residues(row["epitope_project_residues"])
            hits = len(epitope & prior)
            row.update(hotspot_hits_4p5=hits, hotspot_coverage=hits/len(prior),
                       hotspot_fraction_in_interface=hits/len(epitope) if epitope else 0)
        else:
            for field in ("hotspot_hits_4p5", "hotspot_coverage", "hotspot_fraction_in_interface"):
                row[field] = NA


def historical(root):
    records = {}
    comparison = first_existing(root / "results/14_quality/formal_msa_comparison.csv", root / "results/formal_msa_comparison.csv")
    for row in read(comparison):
        condition = "empty" if row["condition"] == "formal" else row["condition"]
        records[(ROUTES[row["route"]], condition)] = dict(row)
    for row in read(root / "docs/route_ensemble_summary.csv"):
        records.setdefault((ROUTES[row["route"]], "empty"), {}).update(row)
    report_path = root / "docs/Phase2_summary.md"
    report = report_path.read_text(encoding="utf-8") if report_path.is_file() else ""
    for match in re.finditer(r"- (native_blind|immunogen_blind)( with formal MSA)?: mean PRODIGY dG ([-0-9.]+) kcal/mol \(SD ([0-9.]+)\)", report):
        route, formal, mean, sd = match.groups()
        record = records.setdefault((ROUTES[route], "formal_msa" if formal else "empty"), {})
        record.update(prodigy_dG_kcal_mol_mean=mean, prodigy_dG_kcal_mol_sd=sd,
                      prodigy_source_path=str(report_path), prodigy_evidence_level="historical_report_only")
    return records


def audit_materials(root, audit, models, summaries, pdbs, prior_path):
    files = [root / p for p in (
        "results/complex_model_summary.csv", "results/interface_contacts.csv",
        "results/model_ranking_formal_msa.csv", "results/prodigy_formal_msa.csv",
        "results/formal_msa_comparison.csv", "docs/route_ensemble_summary.csv",
        "docs/Phase2_summary.md", "MODEL_CARD.md", "REPRODUCIBILITY.md",
        "results.csv", "results_immunogen.csv", "src/rank_formal_models.py",
        "src/rank_formal_msa_models.py", "src/summarize_formal_msa.py",
        "src/analyze_boltz_complexes.py", "src/generate_af3_input.py")]
    for condition in CONDITIONS:
        files.extend(baseline_paths(root, condition).values())
    files.extend((prior_path, audit / "server_raw_model_audit.csv"))
    write(audit / "source_file_metadata.csv", [dict(base("multiple", "multiple", "ALL", p), **{k:v for k,v in metadata(p).items() if k != "source_path"}) for p in sorted(set(files))])
    prior_rows = read(prior_path)
    prior = {int(r["project_residue"]) for r in prior_rows} if prior_rows else None
    if prior is not None and (not prior or min(prior) < 1 or max(prior) > 181):
        raise ValueError("Frozen prior contains invalid project positions")
    prior_record = metadata(prior_path)
    prior_record.update(status="available" if prior else "blocked_missing_frozen_prior",
                        n_residues=len(prior) if prior else NA,
                        note="Read only the frozen CSV; never reconstruct it from source-code ranges or Full181 outputs. No hashes computed.")
    (audit / "frozen_prior_status.json").write_text(json.dumps(prior_record, indent=2), encoding="utf-8")
    contacts = []
    contact_paths = {}
    for condition in CONDITIONS:
        contact_paths[condition] = baseline_paths(root, condition)["contacts"]
        contacts.extend(normalized_legacy(read(contact_paths[condition]), condition, contact_paths[condition]))
    contact_sets = defaultdict(set)
    for row in contacts:
        if int(row["contact_4p5"]):
            contact_sets[key(row)].add(int(row["antigen_project"]))
        if int(row["antigen_project"]) != int(row["antigen_local"]) + 26:
            raise ValueError("Legacy contact coordinate mismatch")
    checks = []
    for row in summaries:
        path = contact_paths[row["condition"]]
        if not path.is_file():
            checks.append(dict(base(row["construct"], row["condition"], row["seed"], path),
                               check="saved_contact_epitope_vs_saved_summary", status="unavailable_contact_table"))
            continue
        expected = residues(row["epitope_project_residues"])
        observed = contact_sets[key(row)]
        good = expected == observed and len(expected) == int(row["epitope_residues_4p5"])
        checks.append(dict(base(row["construct"], row["condition"], row["seed"], path),
                           check="saved_contact_epitope_vs_saved_summary", status="passed" if good else "failed"))
        if not good:
            raise ValueError("Saved baseline contact/summary disagreement")
    write(audit / "actual_data_validation.csv", checks,
          ["construct", "condition", "seed", "source_path", "check", "status"])
    by_condition = Counter(r["condition"] for r in summaries)
    verified_bundles = sum(r["current_raw_prediction_bundle_verified"] for r in models)
    text = (
        "# Baseline audit\n\n"
        "Historical documentation records 15 seeds per construct and condition (60 models). "
        f"The current selected inputs contain {len(summaries)} saved per-model summaries "
        f"(empty={by_condition['empty']}, formal_msa={by_condition['formal_msa']}), "
        f"{len(contacts)} saved contact rows, and {len(pdbs)} PDB files representing "
        f"{len({key(r) for r in pdbs})} unique models. "
        f"The optional server CPU audit verifies {verified_bundles} accessible original PDB/confidence/PAE bundles. "
        "Missing per-model rows remain NA; historical reported completion is separate from current availability.\n\n"
        "Original per-condition table directories take precedence over flattened release fallbacks. "
        "When a flattened summary is used, its formal-MSA identity is checked against saved ranking metrics. "
        "results.csv is only a selected candidate export and never substitutes for an ensemble; "
        "an empty results_immunogen.csv does not imply Immunogen was not run.\n\n"
        "Original ranking: confidence descending; consensus mean pairwise epitope Jaccard descending "
        "then confidence; informed frozen-prior coverage descending then confidence. The older empty-MSA "
        "ranking did not contain consensus ranking. Full181 applies the formal-MSA consensus rule separately "
        "to each condition, with numeric seed ascending as the fixed final tie break. "
        "Prior coverage denominator is the entire frozen prior, not the construct intersection.\n\n"
        "Consensus uses >=50% support: 8/15. Legacy Jaccard assigns 1 to two empty epitopes; "
        "the exports separately mark both_empty and exclude it from nonempty-union summaries. "
        "Incomplete Full181 consensus is provisional only; absent models never count as no-contact models.\n\n"
        "Frozen prior CSV status: " + prior_record["status"] + ". "
        "The historical generator lists ranges, but this is not substituted for the missing frozen file. "
        "Saved old hotspot scores remain labelled historical. No hash was computed.\n\n"
        "The empty-MSA shared 27–160 consensus Jaccard is 1/31=0.032258 from the frozen sets, "
        "consistent with formal_msa_comparison.csv; Phase2_summary.md gives 0.031, which is not the "
        "shared-region recomputation. Geometry uses only actual supplied CPU audit observations; "
        "no subset is extrapolated to 15 models. No old GPU or inference task was launched.\n")
    (audit / "baseline_audit.md").write_text(text, encoding="utf-8")
    return prior


def ensemble(models, history, out):
    groups = defaultdict(list)
    for row in models:
        groups[(row["construct"], row["condition"])].append(row)
    support, pairs, aggregate, consensus_sets = [], [], [], {}
    for (construct, condition), rows in groups.items():
        valid = [r for r in rows if r["observation_available"]]
        sets = {int(r["seed"]): residues(r["epitope_project_residues"]) for r in valid}
        counts = Counter(p for positions in sets.values() for p in positions)
        hist = history.get((construct, condition), {})
        historical_consensus = residues(hist.get("epitope_consensus_ge50pct"))
        if len(valid) == 15:
            consensus = {p for p, n in counts.items() if n >= 8}
            consensus_state = "complete_15_seed_ensemble"
        elif not valid and historical_consensus is not None:
            consensus = historical_consensus
            consensus_state = "historical_summary_only"
        else:
            consensus = None
            consensus_state = "missing_models_official_consensus_unavailable"
        consensus_sets[(construct, condition)] = (consensus, consensus_state)
        start, end = CONSTRUCTS[construct]["project_start"], CONSTRUCTS[construct]["project_end"]
        source = rows[0].get("model_summary_source_path", rows[0]["source_path"])
        evidence = ";".join(sorted({r["evidence_level"] for r in valid or rows}))
        for position in range(1, 182):
            applicable = start <= position <= end
            observed = applicable and bool(valid)
            record = base(construct, condition, "ALL", source)
            record.update(project_position=position, amino_acid=FULL181_SEQUENCE[position-1],
                          local_position=position-start+1 if applicable else NA,
                          residue_in_construct=applicable, n_planned=15,
                          n_completed_or_saved_models=len(valid), n_missing_models=15-len(valid),
                          support_count=counts[position] if observed else NA,
                          support_rate_completed=counts[position]/len(valid) if observed else NA,
                          support_rate_planned_lower_bound=counts[position]/15 if observed else NA,
                          consensus_ge50pct=int(position in consensus) if applicable and consensus is not None else NA,
                          provisional_consensus_completed=int(counts[position] >= len(valid)/2) if observed else NA,
                          consensus_status=consensus_state if applicable else "not_in_construct",
                          evidence_level=evidence,
                          missing_reason="not_in_construct" if not applicable else "" if observed else "No per-model epitope observations available.")
            support.append(record)
        group_pairs = []
        by_seed = {int(r["seed"]): r for r in rows}
        for left, right in itertools.combinations(range(1, 16), 2):
            record = base(construct, condition, f"{left};{right}", source)
            record.update(seed_left=left, seed_right=right,
                          source_path_left=by_seed[left]["source_path"], source_path_right=by_seed[right]["source_path"],
                          **pair_score(sets.get(left), sets.get(right)))
            pairs.append(record)
            group_pairs.append(record)
        row = base(construct, condition, "ALL", source)
        geometry_observations = [r for r in valid if number(r.get('geometry_warning')) is not None]
        row.update(evidence_level=evidence, n_planned=15,
                   n_current_or_saved_per_model_observations=len(valid), n_missing_per_model_observations=15-len(valid),
                   n_historically_reported_models=hist.get("n_models", NA),
                   n_no_contact_models=sum(not x for x in sets.values()) if valid else NA,
                   n_both_empty_jaccard_pairs=sum(p["both_empty"] == 1 for p in group_pairs) if valid else NA,
                   epitope_consensus_ge50pct=joined(consensus), consensus_status=consensus_state,
                   geometry_ensemble_warning_count=sum(number(r['geometry_warning']) for r in geometry_observations) if len(geometry_observations) == 15 else NA,
                   n_geometry_observations=len(geometry_observations),
                   geometry_observed_warning_count=sum(number(r['geometry_warning']) for r in geometry_observations) if geometry_observations else NA,
                   geometry_missing_reason="" if len(geometry_observations) == 15 else "Only available original structures support geometry; never extrapolate to 15 models.")
        prodigy_rows = [r for r in rows if r.get("prodigy_status") != "missing_observation"]
        row["prodigy_evidence_level"] = ";".join(sorted({r["prodigy_evidence_level"] for r in prodigy_rows})) if prodigy_rows else hist.get("prodigy_evidence_level", "missing_observation")
        row["prodigy_source_path"] = ";".join(sorted({r["prodigy_source_path"] for r in prodigy_rows})) if prodigy_rows else hist.get("prodigy_source_path", NA)
        row["prodigy_status_counts"] = json.dumps(dict(Counter(r.get("prodigy_status", "missing_observation") for r in rows)), sort_keys=True)
        row["prodigy_failure_reasons"] = ";".join(sorted({r["prodigy_failure_reason"] for r in rows if r.get("prodigy_failure_reason")})) or NA
        row["prodigy_temperature_C_values"] = joined({number(r.get("prodigy_temperature_C")) for r in rows if number(r.get("prodigy_temperature_C")) is not None}) or NA
        row["prodigy_configured_temperature_C_values"] = joined({number(r.get("prodigy_configured_temperature_C")) for r in rows if number(r.get("prodigy_configured_temperature_C")) is not None}) or NA
        row["geometry_evidence_levels"] = ";".join(sorted({r.get("geometry_evidence_level", NA) for r in geometry_observations})) or NA
        for metric in METRICS:
            metric_stats = stats(r.get(metric) for r in valid)
            for stat, value in metric_stats.items():
                row[f"{metric}_{stat}"] = value
            if not valid and metric+"_mean" in hist:
                row[f"{metric}_mean"] = hist[metric+"_mean"]
                row[f"{metric}_sd"] = hist.get(metric+"_sd", NA)
                row[f"{metric}_n"] = hist.get("n_models", NA)
        for field in ("jaccard_legacy", "jaccard_nonempty_union"):
            for stat, value in stats(p[field] for p in group_pairs).items():
                row[f"pairwise_{field}_{stat}"] = value
        if not valid and hist:
            row["pairwise_jaccard_legacy_mean"] = hist.get("pairwise_epitope_jaccard_mean", NA)
            row["pairwise_jaccard_legacy_sd"] = hist.get("pairwise_epitope_jaccard_sd", NA)
        aggregate.append(row)
    write(out / "three_construct_epitope_support.csv", support)
    write(out / "full181_epitope_support.csv", [r for r in support if r["construct"] == "Full181"])
    write(out / "three_construct_epitope_jaccard.csv", pairs)
    write(out / "full181_epitope_jaccard.csv", [r for r in pairs if r["construct"] == "Full181"])
    write(out / "three_construct_comparison.csv", aggregate)
    write(out / "full181_ensemble_summary.csv", [r for r in aggregate if r["construct"] == "Full181"])
    return groups, consensus_sets


def shared_comparison(groups, consensus_sets, out):
    rows = []
    comparisons = [(a, b, 27, 160) for a,b in itertools.combinations(CONSTRUCTS, 2)]
    comparisons.append(("Full181", "Native155", 27, 181))
    for condition in CONDITIONS:
        for a, b, start, end in comparisons:
            region = set(range(start, end+1))
            aa, bb = groups[(a, condition)], groups[(b, condition)]
            for left, right in itertools.product(aa, bb):
                x = residues(left["epitope_project_residues"]) if left["observation_available"] else None
                y = residues(right["epitope_project_residues"]) if right["observation_available"] else None
                x, y = x & region if x is not None else None, y & region if y is not None else None
                row = base(a+";"+b, condition, f"{left['seed']};{right['seed']}", left["source_path"])
                row.update(comparison_type="all_cross_seed_pairs", construct_left=a, construct_right=b,
                           seed_left=left["seed"], seed_right=right["seed"],
                           source_path_right=right["source_path"], region=f"project{start}-{end}",
                           epitope_left=joined(x), epitope_right=joined(y), **pair_score(x,y))
                rows.append(row)
            x, state_a = consensus_sets[(a, condition)]
            y, state_b = consensus_sets[(b, condition)]
            x, y = x & region if x is not None else None, y & region if y is not None else None
            row = base(a+";"+b, condition, "ALL", aa[0]["source_path"])
            row.update(comparison_type="consensus_ge50pct", construct_left=a, construct_right=b,
                       seed_left="ALL", seed_right="ALL", source_path_right=bb[0]["source_path"],
                       region=f"project{start}-{end}", epitope_left=joined(x), epitope_right=joined(y),
                       consensus_status_left=state_a, consensus_status_right=state_b, **pair_score(x,y))
            rows.append(row)
    write(out / "shared_region_epitope_comparison.csv", rows)


def rank_full(groups, prior, out):
    ranking, selections = [], []
    reps = out / "representatives"
    reps.mkdir(exist_ok=True)
    for condition in CONDITIONS:
        all_rows = groups[("Full181", condition)]
        valid = [r for r in all_rows if r["observation_available"]]
        scores = {}
        for row in valid:
            epitope = residues(row["epitope_project_residues"])
            other = [pair_score(epitope, residues(r["epitope_project_residues"]))["jaccard_legacy"]
                     for r in valid if r["seed"] != row["seed"]]
            scores[int(row["seed"])] = dict(
                confidence_score=number(row.get("confidence_score")),
                consensus_score=statistics.mean(other) if other else None,
                phase1_informed_score=len(epitope & prior)/len(prior) if prior else None)
        order = {}
        for criterion in ("confidence", "consensus", "phase1_informed"):
            field = criterion+"_score"
            order[criterion] = sorted(
                [int(r["seed"]) for r in valid if scores[int(r["seed"])][field] is not None],
                key=lambda s: (-scores[s][field], -(scores[s]["confidence_score"] or 0), s))
        for row in all_rows:
            seed = int(row["seed"])
            record = base("Full181", condition, seed, row["source_path"])
            record.update(status=row["status"], n_completed_models=len(valid), n_planned_models=15,
                          missing_reason=row["missing_reason"], epitope_project_residues=row["epitope_project_residues"],
                          no_contact_model=int(not residues(row["epitope_project_residues"])) if row["observation_available"] else NA)
            iptm = number(row.get("iptm"))
            record["quality_class"] = ("high" if iptm >= .6 else "medium" if iptm >= .45 else "low") if iptm is not None else NA
            for criterion, ordered in order.items():
                record[criterion+"_rank"] = ordered.index(seed)+1 if seed in ordered else NA
                record[criterion+"_score"] = scores.get(seed, {}).get(criterion+"_score", NA)
            record["phase1_informed_status"] = "available" if prior else "blocked_missing_frozen_prior"
            for field in ("prior_source_path", "prior_residue_count", "hotspot_hits_4p5", "hotspot_coverage", "hotspot_fraction_in_interface"):
                record[field] = row.get(field, NA)
            ranking.append(record)
        for criterion, ordered in order.items():
            chosen = next((r for r in valid if int(r["seed"]) == ordered[0]), None) if ordered else None
            row = base("Full181", condition, chosen["seed"] if chosen else NA, chosen["source_path"] if chosen else out / "full181_model_summary.csv")
            row.update(selection=criterion, selection_rule={
                "confidence": "confidence_score descending; numeric seed ascending",
                "consensus": "mean legacy Jaccard to other completed models descending; confidence descending; numeric seed ascending",
                "phase1_informed": "frozen-prior coverage descending; confidence descending; numeric seed ascending"}[criterion],
                selected=bool(chosen), representative_path=NA,
                status="selected" if chosen else "blocked_missing_frozen_prior" if criterion == "phase1_informed" and not prior else "no_rankable_completed_model")
            if chosen:
                candidate = next((Path(chosen[k]) for k in ("pdb_path", "source_pdb", "source_path")
                                  if chosen.get(k) and Path(chosen[k]).is_file() and Path(chosen[k]).suffix.lower() == ".pdb"), None)
                if candidate:
                    destination = reps / f"Full181_{condition}_seed{chosen['seed']}.pdb"
                    shutil.copy2(candidate, destination)
                    row["representative_path"] = str(destination)
                else:
                    row["status"] = "selected_structure_missing"
                row["interface_state"] = "no_contact" if not residues(chosen["epitope_project_residues"]) else "nonempty_epitope"
            selections.append(row)
    write(out / "full181_model_ranking.csv", ranking)
    write(reps / "selection_manifest.csv", selections)
    write(out / "selection_manifest.csv", selections)


def terminal_and_distributions(models, out):
    terminals, distribution, regions = [], [], []
    for row in models:
        available = row["observation_available"]
        epitope = residues(row["epitope_project_residues"]) if available else None
        for metric in METRICS:
            metric_source = row.get("prodigy_source_path", row["source_path"]) if metric.startswith("prodigy_") else row["source_path"]
            record = base(row["construct"], row["condition"], row["seed"], metric_source)
            record.update(metric=metric, value=row.get(metric, NA) if available else NA,
                          evidence_level=row["evidence_level"], missing_reason=row["missing_reason"])
            if metric.startswith("prodigy_"):
                record.update(evidence_level=row.get("prodigy_evidence_level", NA),
                              status=row.get("prodigy_status", NA), missing_reason=row.get("prodigy_failure_reason", NA),
                              temperature_C=row.get("prodigy_temperature_C", NA))
            distribution.append(record)
        for start, end in ((27, 160), (27, 181)):
            applicable = end <= CONSTRUCTS[row["construct"]]["project_end"]
            record = base(row["construct"], row["condition"], row["seed"], row["source_path"])
            record.update(region=f"project{start}-{end}", region_in_construct=applicable,
                          epitope_project_residues=joined(epitope & set(range(start,end+1))) if available and applicable else NA,
                          evidence_level=row["evidence_level"], status="available" if available and applicable else "not_in_construct" if not applicable else "missing_observation")
            regions.append(record)
        if row["construct"] != "Full181":
            continue
        for start, end in ((1,26), (27,160), (161,181)):
            positions = epitope & set(range(start,end+1)) if available else None
            record = base("Full181", row["condition"], row["seed"], row["source_path"])
            record.update(region=f"project{start}-{end}", contact_residue_count=len(positions) if available else NA,
                          fraction_of_full_epitope=len(positions)/len(epitope) if available and epitope else NA,
                          epitope_project_residues=joined(positions),
                          interface_state="missing_observation" if not available else "no_contact_fraction_undefined" if not epitope else "nonempty_epitope",
                          missing_reason=row["missing_reason"])
            terminals.append(record)
    write(out / "terminal_contact_summary.csv", terminals)
    write(out / "three_construct_metric_distribution.csv", distribution)
    write(out / "shared_region_model_epitopes.csv", regions)


def cdr_participation(root, models, out):
    contacts = defaultdict(list)
    for condition in CONDITIONS:
        path = baseline_paths(root, condition)["contacts"]
        for row in normalized_legacy(read(path), condition, path):
            contacts[key(row)].append(row)
    for row in read(out / "full181_interface_contacts.csv"):
        contacts[key(row)].append(row)
    rows = []
    for model in models:
        contact_path = baseline_paths(root, model["condition"])["contacts"] if model["construct"] != "Full181" else out / "full181_interface_contacts.csv"
        valid = model["observation_available"] and contact_path.is_file()
        for chain, region in itertools.product(("H", "L"), ("CDR1", "CDR2", "CDR3", "FR1", "FR2", "FR3", "FR4")):
            selected = {str(r.get("antibody_local", r.get("antibody_local_position")))
                        for r in contacts[key(model)]
                        if str(r.get("contact_4p5", "0")) in ("1", "True", "true")
                        and r.get("antibody_chain") == chain and r.get("antibody_region") == region}
            record = base(model["construct"], model["condition"], model["seed"],
                          contact_path)
            record.update(antibody_chain=chain, antibody_region=region,
                          contacting_residue_count=len(selected) if valid else NA,
                          participates=int(bool(selected)) if valid else NA,
                          evidence_level=model["evidence_level"], missing_reason=model["missing_reason"] if valid or not model["observation_available"] else "Contact table absent; cannot infer CDR participation.")
            rows.append(record)
    write(out / "three_construct_cdr_participation.csv", rows)


def regression(root, groups, consensus_sets, audit):
    checks = []
    source = first_existing(root / "results/14_quality/formal_msa_comparison.csv", root / "results/formal_msa_comparison.csv")
    for saved in read(source):
        condition = "empty" if saved["condition"] == "formal" else saved["condition"]
        construct = ROUTES[saved["route"]]
        rows = [r for r in groups[(construct, condition)] if r["observation_available"]]
        if not rows:
            continue
        sets = [residues(r["epitope_project_residues"]) for r in rows]
        recomputed = {metric+"_"+stat: value for metric in ("iptm", "complex_plddt", "hotspot_coverage", "epitope_residues_4p5")
                      for stat, value in stats(r[metric] for r in rows).items() if stat in ("mean", "sd")}
        jac = stats(pair_score(a,b)["jaccard_legacy"] for a,b in itertools.combinations(sets,2))
        recomputed.update(pairwise_epitope_jaccard_mean=jac["mean"], pairwise_epitope_jaccard_sd=jac["sd"])
        for metric, value in recomputed.items():
            delta = abs(float(saved[metric])-float(value))
            record = base(construct, condition, "ALL", source)
            record.update(metric=metric, saved=saved[metric], recomputed=value,
                          absolute_difference=delta, status="passed" if delta < 1e-12 else "failed")
            checks.append(record)
        actual = consensus_sets[(construct, condition)][0]
        checks.append(dict(base(construct, condition, "ALL", source),
                           metric="epitope_consensus_ge50pct", saved=saved["epitope_consensus_ge50pct"],
                           recomputed=joined(actual), status="passed" if actual == residues(saved["epitope_consensus_ge50pct"]) else "failed"))
    write(audit / "saved_summary_regression.csv", checks,
          list(dict.fromkeys(k for row in checks for k in row)) or
          ["construct", "condition", "seed", "source_path", "metric", "status"])
    if any(r["status"] != "passed" for r in checks):
        raise ValueError("Baseline saved-table regression failed")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--output", type=Path)
    ap.add_argument("--prior", type=Path, help="Existing frozen phase1_prior_deduplicated.csv; never reconstructed")
    args = ap.parse_args()
    root = args.root.resolve()
    out = args.output.resolve() if args.output else root / "results/full181_supplement"
    audit = out / "baseline_audit"
    audit.mkdir(parents=True, exist_ok=True)
    old, saved, pdbs = load_baseline(root, audit)
    prior_path = (args.prior or root / "results/12_complex_prediction/phase1_prior_deduplicated.csv").resolve()
    prior = audit_materials(root, audit, old, saved, pdbs, prior_path)
    models = load_full(out) + old
    add_full_prior(models, prior, prior_path)
    legacy_geometry = keyed_rows(read(out / 'interface_audit/legacy_raw_subset_interface_qc.csv'), "legacy subset geometry")
    for row in models:
        if row['construct'] == 'Full181' or key(row) not in legacy_geometry or row.get('raw_audit_source_path'):
            continue
        qc = legacy_geometry[key(row)]
        add_geometry(row, qc, 'recomputed_existing_original_PDB_subset')
    if len({key(r) for r in models}) != 90:
        raise ValueError("Expected exactly 90 distinct planned construct/condition/seed keys")
    for row in models:
        if row["observation_available"]:
            positions = residues(row["epitope_project_residues"])
            spec = CONSTRUCTS[row["construct"]]
            if positions - set(range(spec["project_start"], spec["project_end"]+1)):
                raise ValueError("Epitope position outside construct")
    write(out / "three_construct_model_comparison.csv", models)
    groups, consensus = ensemble(models, historical(root), out)
    shared_comparison(groups, consensus, out)
    rank_full(groups, prior, out)
    terminal_and_distributions(models, out)
    cdr_participation(root, models, out)
    regression(root, groups, consensus, audit)
    accepted = sum(r["observation_available"] for r in models if r["construct"] == "Full181")
    result = dict(full181_completed=accepted, full181_planned=30,
                  baseline_saved_per_model_observations=sum(r["observation_available"] for r in old),
                  baseline_saved_by_condition={c:sum(r["observation_available"] for r in old if r["condition"] == c) for c in CONDITIONS},
                  baseline_verified_raw_bundles=sum(r["current_raw_prediction_bundle_verified"] for r in old),
                  baseline_prodigy_numeric_observations=sum(number(r.get("prodigy_dG_kcal_mol")) is not None for r in old),
                  baseline_pdb_files=len(pdbs),
                  baseline_unique_pdb_models=len({key(r) for r in pdbs}),
                  frozen_prior_status="available" if prior else "blocked_missing_frozen_prior",
                  note="90 planned keys are not 90 completed models. No hashes or prediction runs performed.")
    (audit / "export_status.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


MODULE_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = MODULE_ROOT / "scripts" / "input_nterm.py"


def load_module():
    assert SCRIPT.is_file(), "input_nterm.py is not implemented"
    spec = importlib.util.spec_from_file_location("v2_input_nterm", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_input_qc_is_exact_181_of_181():
    module = load_module()
    qc = module.validate_inputs(
        MODULE_ROOT / "inputs/reference_181.fasta",
        MODULE_ROOT / "inputs/reference_structure_181.pdb",
    )
    assert qc["status"] == "PASS"
    assert qc["reference_length"] == 181
    assert qc["pdb_chain_a_length"] == 181
    assert qc["sequence_identity"] == "181/181"
    assert qc["continuous_project_numbering"] is True


def test_reference_fasta_uses_fixed_deepsig_project_id():
    header = (MODULE_ROOT / "inputs/reference_181.fasta").read_text(encoding="utf-8").splitlines()[0]
    assert header == ">IL24_project_181"


def test_uniprot_signal_maps_to_project_26_27_with_disclosed_mismatch():
    module = load_module()
    reference = module.read_fasta(MODULE_ROOT / "inputs/reference_181.fasta")
    record = json.loads((MODULE_ROOT / "inputs/uniprot_annotation.json").read_text())
    mapping = module.map_uniprot_signal(reference, record)

    assert mapping["uniprot_accession"] == "Q925S4"
    assert mapping["uniprot_signal_start"] == 1
    assert mapping["uniprot_signal_end"] == 65
    assert mapping["project_offset"] == 39
    assert mapping["project_signal_start"] == 1
    assert mapping["project_signal_end"] == 26
    assert mapping["construct_start"] == 27
    assert mapping["signal_region_exact_match"] is True
    assert mapping["mapped_identity"] == "180/181"
    assert mapping["mismatches"] == [
        {"project_position": 31, "project_aa": "S", "uniprot_position": 70, "uniprot_aa": "F"}
    ]


def test_formal_run_is_module_pinned_and_emits_required_outputs(tmp_path):
    module = load_module()
    result = module.run_analysis(output_root=tmp_path)
    assert result["construct_start"] == 27
    assert (tmp_path / "results/input_qc.json").is_file()
    assert (tmp_path / "results/01_n_terminal/n_terminal_boundary.csv").is_file()
    assert (tmp_path / "figures/n_terminal_hydropathy.png").is_file()
    assert (tmp_path / "reports/step_v2_01_n_terminal.md").is_file()
    assert (tmp_path / "results/01_n_terminal/deepsig_summary.csv").is_file()
    assert (tmp_path / "results/01_n_terminal/deepsig_provenance.json").is_file()
    boundary = (tmp_path / "results/01_n_terminal/n_terminal_boundary.csv").read_text(encoding="utf-8")
    assert "UNAVAILABLE_EXTERNAL_LICENSED_PACKAGE" in boundary
    assert "STRONGLY_CONCORDANT" in boundary


def test_cli_has_no_biological_input_override():
    module = load_module()
    parser = module.build_parser()
    options = {action.dest for action in parser._actions}
    assert "fasta" not in options
    assert "pdb" not in options
    assert "uniprot" not in options


def test_deepsig_prediction_parser_maps_signal_and_chain(tmp_path):
    module = load_module()
    output = tmp_path / "deepsig.gff3"
    output.write_text(
        "IL24_project_181\tDeepSig\tSignal peptide\t1\t26\t0.98\t.\t.\tevidence=ECO:0000256\n"
        "IL24_project_181\tDeepSig\tChain\t27\t181\t.\t.\t.\tevidence=ECO:0000256\n",
        encoding="utf-8",
    )

    prediction = module.parse_deepsig_gff3(output, reference_length=181)

    assert prediction == {
        "id": "IL24_project_181",
        "prediction": "SP",
        "signal_peptide_start": 1,
        "signal_peptide_end": 26,
        "signal_peptide_score": 0.98,
        "cleavage": "26|27",
        "chain_start": 27,
        "chain_end": 181,
    }


def test_deepsig_prediction_parser_rejects_wrong_project_length(tmp_path):
    module = load_module()
    output = tmp_path / "deepsig.gff3"
    output.write_text(
        "IL24_project_181\tDeepSig\tChain\t1\t180\t.\t.\t.\tevidence=ECO:0000256\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="project length"):
        module.parse_deepsig_gff3(output, reference_length=181)


def test_deepsig_success_requires_zero_return_code_and_real_output(tmp_path):
    module = load_module()
    output = tmp_path / "deepsig.gff3"
    output.write_text(
        "IL24_project_181\tDeepSig\tChain\t1\t181\t.\t.\t.\tevidence=ECO:0000256\n",
        encoding="utf-8",
    )
    assert module.deepsig_execution_succeeded(0, output) is True
    assert module.deepsig_execution_succeeded(1, output) is False
    assert module.deepsig_execution_succeeded(0, tmp_path / "missing.gff3") is False


def test_nterminal_evidence_rules():
    module = load_module()
    exact = module.classify_nterminal_evidence(deepsig_cleavage_end=26, uniprot_cleavage_end=26)
    nearby = module.classify_nterminal_evidence(deepsig_cleavage_end=24, uniprot_cleavage_end=26)
    conflict = module.classify_nterminal_evidence(deepsig_cleavage_end=19, uniprot_cleavage_end=26)
    assert exact["evidence_level"] == "STRONGLY_CONCORDANT"
    assert exact["construct_start"] == 27
    assert nearby["evidence_level"] == "CONCORDANT_BOUNDARY_NEIGHBORHOOD"
    assert nearby["construct_start"] == 27
    assert conflict["evidence_level"] == "CONFLICTING"
    assert conflict["construct_start"] == 27

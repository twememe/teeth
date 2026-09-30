#!/usr/bin/env python3
"""Full formal-MSA GPU workflow with strict preflight and isolated outputs."""
from __future__ import annotations
import argparse, csv, hashlib, importlib.metadata, json, os, re, shutil, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
EXPECTED={"torch":"2.6.0+cu124","ImmuneBuilder":"1.2","ANARCI":"2026.2.13.2","boltz":"2.2.1","prodigy-prot":"2.4.0","OpenMM":"8.6.0","pdbfixer":"1.12.0"}
ASSET_ROWS=list(csv.DictReader((ROOT/"models/MODEL_ASSETS.csv").open(encoding="utf-8")))
FORMAL_SCRIPTS=["run_imgt_numbering.py","step11_antibody_structure.py","generate_boltz_inputs.py","run_boltz_formal_msa.sh","analyze_boltz_complexes.py","run_prodigy_formal.py","rank_formal_msa_models.py","summarize_formal_msa.py","build_recomputed_results.py"]
def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""): h.update(block)
    return h.hexdigest()
def resolved_output(value):
    out=(ROOT/value).resolve() if value else (ROOT/"work"/("full_recompute_"+datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"))).resolve()
    allowed=(ROOT/"work").resolve()
    if out==allowed or allowed not in out.parents: raise ValueError("output directory must be a new child of work/")
    return out
def preflight(out):
    problems=[]
    if sys.version_info[:2]!=(3,10): problems.append(f"Python 3.10.x required; found {sys.version.split()[0]}")
    for package,expected in EXPECTED.items():
        try: actual=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError: problems.append(f"missing Python package: {package}=={expected}"); continue
        if actual!=expected: problems.append(f"wrong {package} version: expected {expected}, found {actual}")
    try:
        import torch
        if torch.__version__!="2.6.0+cu124": problems.append(f"wrong torch runtime: {torch.__version__}")
        if torch.version.cuda!="12.4": problems.append(f"wrong torch CUDA runtime: {torch.version.cuda}")
        if not torch.cuda.is_available(): problems.append("CUDA is unavailable; formal Boltz GPU prediction requires CUDA")
    except Exception as exc: problems.append(f"torch import failed: {exc}")
    model_root=Path(os.environ.get("PHASE2_MODEL_CACHE",ROOT/"models")).resolve()
    for row in ASSET_ROWS[:2]:
        path=model_root/row["group"]/row["asset"]
        if not path.is_file(): problems.append(f"missing model asset: {path}")
        elif path.stat().st_size!=int(row["size_bytes"]) or sha256(path)!=row["sha256"]: problems.append(f"invalid model asset size/hash: {path}")
    try:
        import ImmuneBuilder
        trained=Path(ImmuneBuilder.__file__).resolve().parent/"trained_model"
        for row in ASSET_ROWS[2:]:
            path=trained/row["asset"]
            if not path.is_file(): problems.append(f"ImmuneBuilder runtime weight missing: {path}")
            elif path.stat().st_size!=int(row["size_bytes"]) or sha256(path)!=row["sha256"]: problems.append(f"ImmuneBuilder runtime weight invalid: {path}")
    except Exception as exc: problems.append(f"ImmuneBuilder runtime inspection failed: {exc}")
    for rel in ["data/IL24_native_27_181.fasta","data/IL24_immunogen_27_160.fasta","data/VH_VL_clean.fasta","data/phase1_prior_deduplicated.csv","data/msa/msa_manifest.json"]:
        if not (ROOT/rel).is_file(): problems.append(f"missing input: {rel}")
    manifest=ROOT/"data/msa/msa_manifest.json"
    if manifest.is_file():
        try:
            for row in json.loads(manifest.read_text(encoding="utf-8")):
                p=ROOT/row["path"]
                if not p.is_file(): problems.append(f"missing MSA: {row['path']}")
                elif sha256(p)!=row["sha256"]: problems.append(f"MSA hash mismatch: {row['path']}")
        except Exception as exc: problems.append(f"invalid MSA manifest: {exc}")
    hmmer=os.environ.get("PHASE2_HMMER_BIN")
    hmmer_dir=Path(hmmer).resolve() if hmmer else (Path(shutil.which("hmmscan")).resolve().parent if shutil.which("hmmscan") else None)
    if hmmer_dir is None or not (hmmer_dir/"hmmscan").is_file(): problems.append("HMMER hmmscan not found; set PHASE2_HMMER_BIN")
    for name in FORMAL_SCRIPTS:
        if not (ROOT/"src"/name).is_file(): problems.append(f"missing formal step script: src/{name}")
    if out.exists() and any(out.iterdir()): problems.append(f"output directory is not empty: {out}")
    if shutil.disk_usage(ROOT).free<20*1024**3: problems.append("less than documented 20 GiB free disk")
    for p in ROOT.rglob("*"):
        if p.name=="__pycache__" or p.suffix==".pyc" or p.name.endswith(".partial") or p.name.endswith(".lock"): problems.append(f"temporary artifact present: {p.relative_to(ROOT)}")
    scan=[ROOT/"README.md",ROOT/"docs",ROOT/"models",ROOT/"scripts",ROOT/"src",ROOT/"run.sh",ROOT/"run_full.sh"]
    absolute=re.compile(r"/(?:home|mnt)/|[A-Za-z]:\\\\")
    for base in scan:
        files=[base] if base.is_file() else list(base.rglob("*"))
        for p in files:
            if p.is_file() and absolute.search(p.read_text(errors="ignore")): problems.append(f"absolute path dependency found: {p.relative_to(ROOT)}")
    return sorted(set(problems))
def commands():
    py=sys.executable
    return [[py,"-B",str(ROOT/"src/run_imgt_numbering.py")],[py,"-B",str(ROOT/"src/step11_antibody_structure.py")],[py,"-B",str(ROOT/"src/generate_boltz_inputs.py")],["bash",str(ROOT/"src/run_boltz_formal_msa.sh")],[py,"-B",str(ROOT/"src/analyze_boltz_complexes.py"),"--set","formal_msa"],[py,"-B",str(ROOT/"src/run_prodigy_formal.py"),"--set","formal_msa"],[py,"-B",str(ROOT/"src/rank_formal_msa_models.py")],[py,"-B",str(ROOT/"src/summarize_formal_msa.py")],[py,"-B",str(ROOT/"src/build_recomputed_results.py")]]
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--preflight-only",action="store_true"); ap.add_argument("--output-dir"); ap.add_argument("--self-test-failure",action="store_true",help=argparse.SUPPRESS); a=ap.parse_args()
    if a.self_test_failure: return subprocess.run([sys.executable,"-c","raise SystemExit(23)"],check=False).returncode
    try: out=resolved_output(a.output_dir)
    except ValueError as exc: print(f"PREFLIGHT_FAILED\n- {exc}"); return 2
    problems=preflight(out)
    print(f"OUTPUT_TARGET={out.relative_to(ROOT)}")
    if problems: print("PREFLIGHT_FAILED\n"+"\n".join("- "+p for p in problems)); return 1
    print(f"PREFLIGHT_OK output={out.relative_to(ROOT)}")
    if a.preflight_only: return 0
    out.mkdir(parents=True,exist_ok=False)
    env=os.environ.copy(); env["PHASE2_OUTPUT_ROOT"]=str(out); env["PYTHONDONTWRITEBYTECODE"]="1"
    if not env.get("PHASE2_HMMER_BIN") and shutil.which("hmmscan"):
        env["PHASE2_HMMER_BIN"]=str(Path(shutil.which("hmmscan")).resolve().parent)
    with (out/"run.log").open("w",encoding="utf-8") as lf:
        lf.write("models=ImmuneBuilder 1.2; Boltz-1 2.2.1; PRODIGY 2.4.0\n")
        for cmd in commands():
            lf.write(f"START {datetime.now(timezone.utc).isoformat()} cmd={cmd!r}\n"); lf.flush()
            result=subprocess.run(cmd,cwd=ROOT,env=env,stdout=lf,stderr=subprocess.STDOUT,check=False)
            lf.write(f"END {datetime.now(timezone.utc).isoformat()} exit={result.returncode}\n"); lf.flush()
            if result.returncode: return result.returncode
    return 0
if __name__=="__main__": raise SystemExit(main())

#!/usr/bin/env python3
import hashlib, json, platform, subprocess
from pathlib import Path
import torch
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"/"17_delivery"; OUT.mkdir(parents=True,exist_ok=True)
env={"python":platform.python_version(),"platform":platform.platform(),"torch":torch.__version__,"cuda_runtime":torch.version.cuda,"cuda_available":torch.cuda.is_available(),"gpu":torch.cuda.get_device_name(0) if torch.cuda.is_available() else None}
(OUT/"environment_final.json").write_text(json.dumps(env,indent=2)+"\n")
(OUT/"pip_freeze_final.txt").write_text(subprocess.check_output([str(ROOT/"venv/bin/pip"),"freeze"],text=True))
targets=[ROOT/"data/phase2_input_manifest_sha256.csv",ROOT/"models/boltz/boltz1_conf.ckpt",ROOT/"models/boltz/ccd.pkl",ROOT/"results/13_interface_analysis/formal/complex_model_summary.csv",ROOT/"results/15_affinity/prodigy_formal.csv",ROOT/"docs/Phase2_summary.md"]
lines=["file,bytes,sha256"]
for p in targets:
    h=hashlib.sha256(p.read_bytes()).hexdigest(); lines.append(f"{p.relative_to(ROOT)},{p.stat().st_size},{h}")
(OUT/"final_artifact_manifest_sha256.csv").write_text("\n".join(lines)+"\n")
print(json.dumps(env))

#!/usr/bin/env python3
"""Download, install and verify the exact large model assets."""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, os, shutil, sys
from pathlib import Path
from urllib.request import Request, urlopen
ROOT=Path(__file__).resolve().parents[1]
ASSETS=list(csv.DictReader((ROOT/"models/MODEL_ASSETS.csv").open(encoding="utf-8")))
def digest(path):
    h=hashlib.sha256(); n=0
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""): h.update(block); n+=len(block)
    return n,h.hexdigest()
def valid(path,row):
    return path.is_file() and digest(path)==(int(row["size_bytes"]),row["sha256"])
def immune_runtime_dir():
    spec=importlib.util.find_spec("ImmuneBuilder")
    return Path(spec.origin).resolve().parent/"trained_model" if spec and spec.origin else None
def runtime_path(row,cache):
    if row["group"]=="boltz": return cache/"boltz"/row["asset"]
    target=immune_runtime_dir(); return target/row["asset"] if target else None
def download(row,destination):
    destination.parent.mkdir(parents=True,exist_ok=True); partial=destination.with_name(destination.name+".partial"); partial.unlink(missing_ok=True)
    urls=[row["primary_url"]]+([row["backup_url"]] if row["backup_url"] else [])
    for url in urls:
        try:
            with urlopen(Request(url,headers={"User-Agent":"phase2-release/2.0"}),timeout=120) as src, partial.open("wb") as dst:
                while block:=src.read(1024*1024): dst.write(block)
            if valid(partial,row): os.replace(partial,destination); return
            print(f"size/hash mismatch from {url}",file=sys.stderr)
        except Exception as exc: print(f"download failed {url}: {exc}",file=sys.stderr)
        partial.unlink(missing_ok=True)
    raise RuntimeError(f"unable to obtain verified asset: {row['asset']}")
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--cache-dir",type=Path,default=ROOT/"models"); ap.add_argument("--verify-only",action="store_true"); ap.add_argument("--list",action="store_true"); ap.add_argument("--install-immunebuilder",action="store_true"); a=ap.parse_args(); cache=a.cache_dir.resolve()
    if a.list:
        for row in ASSETS:
            print(f"{row['asset']}\t{row['size_bytes']}\t{row['sha256']}\t{row['primary_url']}\tcache={cache/row['group']/row['asset']}\truntime={runtime_path(row,cache)}")
        return 0
    failures=[]
    if a.verify_only:
        for row in ASSETS:
            path=runtime_path(row,cache)
            if path is None: failures.append(f"missing runtime package for asset: {row['asset']} (ImmuneBuilder not installed)")
            elif not valid(path,row): failures.append(f"missing or invalid runtime asset: {path}")
            else: print(f"OK runtime={path}")
    else:
        for row in ASSETS:
            cached=cache/row["group"]/row["asset"]
            try:
                if not valid(cached,row): download(row,cached)
                print(f"OK cache={cached}")
            except Exception as exc: failures.append(str(exc))
        if a.install_immunebuilder and not failures:
            target=immune_runtime_dir()
            if target is None: failures.append("ImmuneBuilder is not installed; cannot install runtime weights")
            else:
                target.mkdir(parents=True,exist_ok=True)
                for row in [x for x in ASSETS if x["group"]=="immunebuilder"]:
                    source=cache/row["group"]/row["asset"]; partial=target/(row["asset"]+".partial")
                    shutil.copyfile(source,partial)
                    if not valid(partial,row): partial.unlink(missing_ok=True); failures.append(f"install verification failed: {row['asset']}"); continue
                    os.replace(partial,target/row["asset"]); print(f"INSTALLED runtime={target/row['asset']}")
    if failures: print("\n".join(failures),file=sys.stderr); return 1
    return 0
if __name__=="__main__": raise SystemExit(main())

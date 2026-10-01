from __future__ import annotations

import json
import platform
from importlib import metadata

import torch


PACKAGES = [
    "discotope3",
    "torch",
    "torch-geometric",
    "biopython",
    "biotite",
    "biotraj",
    "xgboost",
    "pygam",
    "numpy",
    "pandas",
    "joblib",
    "scikit-learn",
]

payload = {
    "python": platform.python_version(),
    "packages": {name: metadata.version(name) for name in PACKAGES},
    "torch_cuda_build": torch.version.cuda,
    "cuda_available": torch.cuda.is_available(),
    "cuda_device_count": torch.cuda.device_count(),
    "cuda_device_name": (
        torch.cuda.get_device_name(0) if torch.cuda.is_available() else None
    ),
}
print(json.dumps(payload, indent=2, sort_keys=True))

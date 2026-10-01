from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from discotope3 import esm
from discotope3 import main as discotope_main


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--models-dir", type=Path, required=True)
    args = parser.parse_args()

    models_dir = args.models_dir.resolve()
    # The official CLI initializes this module-global logger in main().
    # Reproduce that context for a direct loader smoke test without patching upstream.
    discotope_main.log = logging.getLogger("discotope3.main")
    xgb_models = discotope_main.load_models(str(models_dir), num_models=100)
    gam_len = discotope_main.load_gam_model(
        str(models_dir / "gam_len_to_mean.pkl")
    )
    gam_surface = discotope_main.load_gam_model(
        str(models_dir / "gam_surface_to_std.pkl")
    )
    esm_model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()

    payload = {
        "models_dir": str(models_dir),
        "xgboost_models_loaded": len(xgb_models),
        "gam_len_model_class": type(gam_len).__name__,
        "gam_surface_model_class": type(gam_surface).__name__,
        "esm_if1_model_class": type(esm_model).__name__,
        "esm_if1_parameter_count": sum(p.numel() for p in esm_model.parameters()),
        "esm_alphabet_size": len(alphabet),
        "status": "PASS",
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import json
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from euv_playground.multilayer import reflectivity_vs_wavelength
from euv_playground.resist import shot_noise

out = ROOT / "data"; out.mkdir(exist_ok=True)
lam = np.linspace(12.5, 14.5, 401)
(out/"reflectivity.json").write_text(json.dumps({"wavelength_nm": lam.tolist(), "reflectivity": reflectivity_vs_wavelength(lam).tolist()}, indent=2)+"\n")
doses = np.linspace(5, 100, 20); sizes = np.array([5, 7, 10, 15, 20])
_, photons, noise = shot_noise(doses[:, None], sizes[None, :])
(out/"shot_noise.json").write_text(json.dumps({"dose_mJ_cm2": doses.tolist(), "feature_size_nm": sizes.tolist(), "photons_per_pixel": photons.tolist(), "noise_fraction": noise.tolist()}, indent=2)+"\n")
print(f"Wrote data to {out}")

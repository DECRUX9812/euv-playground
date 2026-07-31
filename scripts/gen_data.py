#!/usr/bin/env python3
import json
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from euv_playground.multilayer import reflectivity_vs_wavelength
from euv_playground.resist import shot_noise
from euv_playground.imaging import (abbe_image_2d, contact_holes,
                                    contrast_vs_pitch, dense_lines,
                                    isolated_line)

out = ROOT / "data"; out.mkdir(exist_ok=True)
lam = np.linspace(12.5, 14.5, 401)
(out/"reflectivity.json").write_text(json.dumps({"wavelength_nm": lam.tolist(), "reflectivity": reflectivity_vs_wavelength(lam).tolist()}, indent=2)+"\n")
doses = np.linspace(5, 100, 20); sizes = np.array([5, 7, 10, 15, 20])
_, photons, noise = shot_noise(doses[:, None], sizes[None, :])
(out/"shot_noise.json").write_text(json.dumps({"dose_mJ_cm2": doses.tolist(), "feature_size_nm": sizes.tolist(), "photons_per_pixel": photons.tolist(), "noise_fraction": noise.tolist()}, indent=2)+"\n")

n, dx = 128, 2.0
masks = {
    "dense_lines": dense_lines(n, dx, 64),
    "isolated_line": isolated_line(n, dx, 32),
    "contact_holes": contact_holes(n, dx, 64, 16),
}
cases = {
    "ideal": (), "defocus": ((2, 0, 0.8),),
    "coma": ((3, 1, 0.8),), "astigmatism": ((2, 2, 0.8),),
}
images = {name: {case: abbe_image_2d(mask, dx, zernike_terms=terms).tolist()
                 for case, terms in cases.items()} for name, mask in masks.items()}
(out/"imaging_2d.json").write_text(json.dumps(
    {"dx_nm": dx, "shape": [n, n], "images": images}, separators=(",", ":"))+"\n")

pitches = np.arange(14.0, 81.0, 2.0)
curves = {"sigma_0": contrast_vs_pitch(256, 1.0, pitches, 0.0),
          "sigma_1": contrast_vs_pitch(256, 1.0, pitches, 1.0)}
(out/"resolution_curve.json").write_text(json.dumps(
    {"wavelength_nm": 13.5, "na": 0.33, "curves": curves}, indent=2)+"\n")
print(f"Wrote data to {out}")

#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from euv_playground.multilayer import reflectivity_vs_wavelength
from euv_playground.resist import shot_noise

lam = np.linspace(12.5, 14.5, 4001); r = reflectivity_vs_wavelength(lam)
i = int(np.argmax(r)); support = np.flatnonzero(r >= r[i]/2)
fwhm = lam[support[-1]]-lam[support[0]]
_, photons, noise = shot_noise(30, 10)
checks = [
    ("peak reflectivity", .65 <= r[i] <= .75, f"{100*r[i]:.4f}%"),
    ("peak wavelength", abs(lam[i]-13.5) <= .3, f"{lam[i]:.4f} nm"),
    ("reflectivity FWHM", .3 <= fwhm <= .8, f"{fwhm:.4f} nm"),
    ("photons/10 nm pixel", 15 <= photons <= 25, f"{photons:.4f}"),
    ("shot noise", .18 <= noise <= .27, f"{100*noise:.4f}%"),
]
for name, ok, value in checks:
    print(f"{'PASS' if ok else 'FAIL'}: {name}: {value}")
raise SystemExit(0 if all(x[1] for x in checks) else 1)

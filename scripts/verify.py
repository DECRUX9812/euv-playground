#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from euv_playground.multilayer import reflectivity_vs_wavelength
from euv_playground.resist import shot_noise
from euv_playground.imaging import (abbe_image_2d, contrast,
                                    contrast_vs_defocus, contrast_vs_pitch,
                                    dense_lines)

lam = np.linspace(12.5, 14.5, 4001); r = reflectivity_vs_wavelength(lam)
i = int(np.argmax(r)); support = np.flatnonzero(r >= r[i]/2)
fwhm = lam[support[-1]]-lam[support[0]]
_, photons, noise = shot_noise(30, 10)

def measured_cutoff(sigma):
    pitches = np.arange(14.0, 51.0, 0.5)
    curve = contrast_vs_pitch(256, 1.0, pitches, sigma)
    for (p0, c0), (p1, c1) in zip(curve, curve[1:]):
        if c0 < .1 <= c1:
            return p0 + (.1-c0)*(p1-p0)/(c1-c0)
    return float("nan")

coherent_cutoff = measured_cutoff(0.0)
partial_cutoff = measured_cutoff(1.0)
defocus_values = [c for _, c in contrast_vs_defocus(256, 1.0, 70)]
period = 64
period_image = abbe_image_2d(dense_lines(256, 1.0, period), 1.0)
period_error = float(np.max(np.abs(period_image - np.roll(period_image, period, axis=1))))
checks = [
    ("peak reflectivity", .65 <= r[i] <= .75, f"{100*r[i]:.4f}%"),
    ("peak wavelength", abs(lam[i]-13.5) <= .3, f"{lam[i]:.4f} nm"),
    ("reflectivity FWHM", .3 <= fwhm <= .8, f"{fwhm:.4f} nm"),
    ("photons/10 nm pixel", 15 <= photons <= 25, f"{photons:.4f}"),
    ("shot noise", .18 <= noise <= .27, f"{100*noise:.4f}%"),
    ("coherent imaging cutoff", 34.7 <= coherent_cutoff <= 47.0,
     f"{coherent_cutoff:.2f} nm"),
    ("sigma=1 imaging cutoff", 17.4 <= partial_cutoff <= 23.6,
     f"{partial_cutoff:.2f} nm"),
    ("defocus contrast monotonic", defocus_values[0] > defocus_values[1] > defocus_values[2],
     " > ".join(f"{value:.4f}" for value in defocus_values)),
    ("imaging period preservation", period_error < 2e-3,
     f"max period-shift error {period_error:.3g}"),
]
for name, ok, value in checks:
    print(f"{'PASS' if ok else 'FAIL'}: {name}: {value}")
raise SystemExit(0 if all(x[1] for x in checks) else 1)

# Writeup 01 — Mo/Si multilayer: why EUV is 13.5 nm

**Date:** 2026-07-31
**Tool:** `euv-playground` — Parratt recursion (exact Fresnel stack), 3 independent implementations
**Stack:** N=50 bilayer pairs · d=6.90 nm · Γ = d_Mo/d = 0.40 · σ = 0.3 nm roughness · θ = 0° (normal incidence)

## What we computed

| Quantity | Textbook expectation | **Computed** |
|---|---|---|
| Peak reflectivity | ~70% (measured CXRO values) | **72.83%** |
| Peak wavelength | 13.5 nm (the EUV standard) | **13.45 nm** |
| R @ 13.5 nm | — | 72.47% |
| FWHM (R vs λ) | ~0.5 nm | **0.65 nm** |
| R with zero roughness (σ=0) | — | 73.62% |
| R @ 6° incidence | — | 69.18% |
| R with N=20 pairs | — | 57.99% |
| Shot noise @ 30 mJ/cm², 10 nm pixel | ~20 photons → ±22% | **20.4 photons → ±22.1%** |

## What the numbers tell us

- **The Bragg condition picks the wavelength.** With d = 6.9 nm, the peak lands at
  13.45 nm — 2d·sinθ ≈ λ at normal incidence. The entire EUV industry runs at
  13.5 nm because that's where Mo/Si multilayers reflect best: Mo and Si have
  their absorption edges bracketing this wavelength, so the stack gets maximum
  index contrast with manageable absorption. It's a materials coincidence the
  whole industry is built on.
- **Γ = 0.4 is a compromise.** Mo gives the index contrast that makes the stack
  reflect, but Mo is *absorbing* at 13.5 nm. Make the Mo layers thicker and the
  mirror eats its own light; make them thinner and the contrast dies. 0.4 is the
  sweet spot — thin enough to reflect, thick enough to stay invisible.
- **Roughness is the enemy.** Interface roughness σ = 0.3 nm (a real-world value
  for state-of-the-art deposition) costs ~0.8% peak reflectivity (73.62% → 72.83%).
  In a 10-mirror optical train that compounds to a ~7% throughput loss — before
  you've lost the other 99% to absorption.
- **50 pairs is enough.** N=20 gives 58%, N=50 gives 72.8% — the gain from more
  pairs saturates because the bottom layers barely see light. This is why real
  EUV mirrors are ~40–60 pairs, not thousands.
- **Off-resonance, the mirror is a window.** Move 0.5 nm off peak and R collapses
  from 72% to 23%. EUV optics have zero tolerance for wavelength drift — the
  tin-plasma source must be stabilized to a fraction of a percent.

## The shot-noise punchline

Each 13.5 nm photon carries 92 eV. A typical 30 mJ/cm² resist dose delivers only
**20 photons per 10 nm pixel** — Poisson statistics then give ±22% intensity
noise. About **one in five of the "ink" that writes a modern transistor is
random chance.** Shrink the feature and it gets worse. This is not a
manufacturing nuisance; it's a fundamental counting-statistics wall, and it's
the reason EUV scaling conversations always end in stochastics.

## Verification

- Physics core verified by `pytest` (7/7) and `scripts/verify.py` (5/5 PASS)
  against published CXRO/LBNL reference values.
- Three independent implementations (Python core, JS prototype, in-page JS) agree
  on peak reflectivity to 4 significant figures: **72.83% @ 13.45 nm**.
- Optical constants: Mo δ=0.0777 β=0.00635, Si δ=0.00255 β=0.00211 (13.5 nm,
  CXRO/LBNL tables).

## Raw data

- `data/reflectivity.json` (R vs λ grid) · `data/shot_noise.json` (dose × feature grid)
- `plots/reflectivity_wavelength.png` · `plots/resist_and_shot_noise.png`

## Run it yourself

```bash
git clone https://github.com/DECRUX9812/euv-playground && cd euv-playground
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/pytest -q && .venv/bin/python scripts/verify.py
# or just open web/index.html in a browser — no server, no deps
```

## Next

- 02 — Abbe-method imaging: mask → pupil (Zernike aberrations) → image
- 03 — Aberration budget: how defocus/coma/astigmatism reshape a printed line

# EUV Playground — open-source EUV lithography simulator

Interactive, defensible-physics simulation of extreme ultraviolet (13.5 nm)
lithography. Built from real equations, verified against published reference
data. No NDAs, no ASML — just optics you can compute.

## Why EUV is weird (30 seconds)

Chips are printed by shining light through a mask onto photoresist. Shorter
wavelength → smaller features. The industry stalled at 193 nm and the only way
forward is **EUV at 13.5 nm**. At that wavelength:

- **Everything absorbs** — no lenses exist, only mirrors, and the whole tool
  runs in vacuum.
- **Mirrors are multilayer stacks** — 40–60 pairs of alternating Mo/Si layers,
  ~7 nm per pair, engineered to Bragg-reflect at exactly 13.5 nm. Peak
  reflectivity is only ~70% (vs 99.9% for visible-light mirrors), so a scanner's
  ~10 mirrors throw away ~99% of the source light.
- **Shot noise is brutal** — each 13.5 nm photon carries 92 eV. A 30 mJ/cm²
  resist dose is only ~20 photons per 10 nm feature; Poisson noise on 20
  photons is ±22%. Roughly a fifth of the "ink" that writes a modern transistor
  is random chance.

## Modules

| # | Module | Physics | Verification target |
|---|--------|---------|---------------------|
| 1 | `multilayer.py` | Mo/Si multilayer reflectivity via Parratt recursion (Fresnel stack), R vs wavelength & incidence angle | Peak R ≈ 70% at 13.5 nm, normal incidence, for N=50, d≈6.9 nm, Γ≈0.4 (CXRO/LBNL published values) |
| 2 | `zernike.py` | Zernike polynomial basis + wavefront maps (defocus, coma, astigmatism, spherical) | Orthonormal on unit disk; standard OSA/ANSI indexing |
| 3 | `imaging.py` | Scalar diffraction imaging: mask → pupil (with aberrations) → image, Abbe method with partial coherence | Coherent limit reproduces classical diffraction result; symmetric patterns stay symmetric |
| 4 | `resist.py` | Resist contrast curve (dose vs remaining thickness), dose-to-clear, photon shot noise → line-edge roughness estimate | ~20 photons / 10 nm pixel at 30 mJ/cm² (Poisson: σ/μ ≈ 22%) |

## Reference numbers (from published CXRO / literature data)

- Mo/Si bilayer period d ≈ 6.9–7.0 nm for 13.5 nm normal incidence
  (Bragg condition 2d·sinθ = m·λ with θ ≈ 90° ⇒ d ≈ λ/2 ≈ 6.75 nm).
- Thickness ratio Γ = d_Mo / d ≈ 0.4.
- Peak normal-incidence reflectivity at 13.5 nm: ~69–75% (theoretical),
  ~69–70% (measured). FWHM of R vs λ: ~0.5 nm.
- Photon energy at 13.5 nm: E = 1240 eV·nm / 13.5 nm ≈ 91.8 eV.
- Photon density at dose D: photons/cm² = D / (91.8 eV × 1.602e-19 J/eV).

## Layout

```
euv_playground/        # Python package (physics core)
  multilayer.py        # Parratt recursion reflectivity
  zernike.py           # Zernike polynomials
  imaging.py           # Abbe-method scalar imaging
  resist.py            # contrast curves + shot noise
tests/                 # pytest suite, physics verified against reference numbers
web/                   # self-contained interactive HTML artifact (pure JS)
plots/                 # generated reference plots (PNG/SVG)
data/                  # generated data JSON (for the web artifact)
scripts/               # CLI entrypoints: gen_plots.py, gen_data.py, verify.py
README.md
```

## Commands

```bash
.venv/bin/pip install -r requirements.txt
.venv/bin/pytest -q
.venv/bin/python scripts/gen_plots.py    # writes plots/*.png
.venv/bin/python scripts/gen_data.py     # writes data/*.json for web artifact
.venv/bin/python scripts/verify.py       # prints pass/fail vs reference numbers
```

Physics core is `numpy` + `scipy` only. Plots via `matplotlib`. No
AI-slop: every curve must match a published reference value within tolerance.

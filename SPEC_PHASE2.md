# SPEC — Phase 2: Interactive Imaging with Zernike Aberrations

Extends the EUV Playground with mask → pupil → image simulation (Abbe method),
2D masks, Zernike aberration control, and a live interactive imaging demo.

## Physics contract (authoritative — verify against these, do not weaken)

All simulations: wavelength λ = 13.5 nm, NA = 0.33 (current EUV scanners are
0.33–0.55; we use 0.33).

- **Resolution limit (dense L/S, duty 0.5):** for a binary grating only odd
  orders (±1/p, ±3/p, …) exist; the image needs the ±1 orders to pass the
  pupil. Coherent (σ=0): cutoff pitch p_c ≈ λ/NA ≈ **40.9 nm** — below this,
  image contrast → 0. Partially coherent (Abbe sum over source tilts to σ·NA):
  p_c ≈ λ/(NA·(1+σ)) → with σ=1, p_c ≈ **20.5 nm**. Verify contrast(pitch)
  drops below 0.1 at approximately these cutoffs (±15%).
- **Defocus:** contrast decreases monotonically with |Z20 coefficient| for a
  fixed pitch above cutoff. Defocus does NOT shift the pattern position.
- **Astigmatism:** Z22 (axis 0°) degrades vertical lines more than horizontal;
  Z2-2 (axis 45°) breaks symmetry diagonally. Rotating the term's axis rotates
  the degradation.
- **Aberration-free dense L/S:** image intensity is periodic with the mask
  period (magnification 1) and contrast > 0.7 for pitch ≥ 2× cutoff.
- **Contact-hole array (2D):** ideal case symmetric; coma (Z31) tilts/blurs
  asymmetrically; spherical (Z40) causes isotropic blur + sidelobes.

## File ownership

- **Codex agent:** `euv_playground/`, `tests/`, `scripts/`, `data/`, `plots/`
- **Web agent:** `web/index.html` ONLY (extend in place — do not touch Python)

## Python deliverables (Codex)

1. Extend `euv_playground/imaging.py`:
   - `abbe_image_2d(mask2d, dx, wavelength=13.5, numerical_aperture=0.33,
     zernike_terms=(), sigma_c=0.0, source_points=9)` — 2D FFT-based Abbe
     imaging, circular pupil (r ≤ 1), phase from existing `phase_map` on 2D
     polar grids. Partial coherence: average over a small 2D grid of source
     tilts within σ_c·NA (use a rotated/uniform sampling, not a full square
     grid, for symmetry).
   - 2D mask builders: `dense_lines(n, dx, pitch, duty=0.5)` (vertical),
     `isolated_line(n, dx, width)`, `contact_holes(n, dx, pitch, hole_r)`.
   - `contrast(intensity)` = (max−min)/(max+min) over the image region.
   - `contrast_vs_pitch(n, dx, pitches, sigma_c, zernike_terms=())` →
     list of (pitch, contrast).
   - `contrast_vs_defocus(...)` — contrast for fixed pitch as |Z20| grows.
2. `tests/test_imaging2d.py` asserting the physics contract above:
   - coherent cutoff: contrast < 0.1 at p = 30 nm; contrast > 0.5 at p = 70 nm
   - σ=1 cutoff: contrast < 0.1 at p = 16 nm; contrast > 0.3 at p = 26 nm
   - monotonic contrast decrease with |defocus| (3 points, e.g. Z20 = 0, 0.25,
     0.5 waves: c0 > c1 > c2)
   - aberration-free period preservation (peak spacing ≈ pitch)
   - astigmatism axis rotation: contrast(vertical lines, Z22) <
     contrast(horizontal lines, Z22) for same coefficient
3. `scripts/gen_plots.py` — add:
   - `imaging_2d_grid.png`: 2×2 panel of dense L/S images (ideal, defocus,
     coma, astigmatism) as 2D intensity maps
   - `resolution_curve.png`: contrast vs pitch for σ=0 and σ=1, cutoff markers
   - `pupil_phases.png`: pupil phase maps for defocus/astigmatism/coma/spherical
4. `scripts/gen_data.py` — add `data/imaging_2d.json` (image grids for the
   three masks at ideal + defocus + coma + astigmatism, downsampled to ≤128²)
   and `data/resolution_curve.json` (pitch/contrast for σ=0 and σ=1).
5. `scripts/verify.py` — add PASS/FAIL entries: coherent cutoff in
   [34.7, 47.0] nm; σ=1 cutoff in [17.4, 23.6] nm; defocus monotonic; period
   preservation. Keep all existing checks green.

## Web deliverables (opencode-go)

Add **DEMO 03 · IMAGING** to `web/index.html` (after Demo 2, before the
"Why EUV is weird" section), same dark theme, no external deps:

- Mask selector (radio buttons): dense L/S, isolated line, contact holes.
- Zernike sliders: defocus, astigmatism (axis 0°), coma, spherical — each
  −0.5…+0.5 waves, default 0. Plus a partial-coherence σ slider 0–1 (default 0).
- Two canvases side by side: **pupil phase map** (color wheel: hue = phase,
  brightness = |r|≤1) and **image intensity** (grayscale heatmap).
- Pure JS FFT: implement a radix-2 complex FFT + separable 2D transform
  (row FFTs then column FFTs). The existing complex helpers (cadd/csub/cmul…)
  in the file can be extended — keep them.
- Readout: image contrast (computed from the JS image), and the pitch being
  simulated (e.g. "pitch 60 nm — 1.47× cutoff").
- Defaults (dense L/S, all aberrations 0, σ=0, pitch ~60 nm) must show a clean
  periodic image with contrast ≥ 0.5. Moving defocus to 0.5 waves must visibly
  blur it. Below-coherent-cutoff pitch (~35 nm) must visibly wash out.
- 2D FFT is the classic failure point: verify row/column separability and
  that an aberration-free grating image has exactly the mask period.
- Keep DEMO 1 and DEMO 2 byte-for-byte functional (do not break existing
  sliders/plots).

## Acceptance

`pytest -q` all green (old + new), `scripts/verify.py` all PASS, web page
loads with zero console errors, Demo 3 sliders live-update both canvases,
numbers cross-check between Python core and in-page JS where comparable.

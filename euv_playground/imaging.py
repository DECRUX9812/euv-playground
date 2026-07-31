"""Scalar Abbe imaging of one- and two-dimensional EUV masks."""

from __future__ import annotations
import numpy as np
from .zernike import phase_map


def line_space_grating(x, period, duty_cycle=0.5):
    """Binary amplitude grating centered symmetrically about x=0."""
    x = np.asarray(x, float)
    return (np.abs(((x + period/2) % period)-period/2) < duty_cycle*period/2).astype(float)


def abbe_image(mask, dx, wavelength=13.5, numerical_aperture=0.33,
               zernike_terms=(), sigma_c=0.0, source_points=21):
    """Return image-plane intensity from a sampled 1-D amplitude mask.

    ``dx`` and wavelength share units. Partial coherence is modeled as an
    incoherent average over uniformly sampled source tilts across sigma_c NA.
    """
    mask = np.asarray(mask, complex)
    if mask.ndim != 1 or not 0 <= sigma_c <= 1:
        raise ValueError("mask must be 1-D and sigma_c must lie in [0, 1]")
    n = mask.size
    freq = np.fft.fftshift(np.fft.fftfreq(n, dx))
    cutoff = numerical_aperture/wavelength
    pupil_x = freq/cutoff
    pupil = np.abs(pupil_x) <= 1
    phase = phase_map(np.abs(pupil_x), np.where(pupil_x >= 0, 0.0, np.pi), zernike_terms)
    phase_factor = pupil*np.exp(1j*phase)
    x = (np.arange(n)-n//2)*dx
    tilts = np.array([0.0]) if sigma_c == 0 else np.linspace(-sigma_c*cutoff, sigma_c*cutoff, source_points)
    intensity = np.zeros(n)
    for tilt in tilts:
        illuminated = mask*np.exp(2j*np.pi*tilt*x)
        spectrum = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(illuminated)))
        field = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(spectrum*phase_factor)))
        intensity += np.abs(field)**2
    return intensity/len(tilts)


def image_grating(x, period, duty_cycle=0.5, **kwargs):
    """Convenience wrapper returning a grating mask and its image intensity."""
    mask = line_space_grating(x, period, duty_cycle)
    return mask, abbe_image(mask, float(x[1]-x[0]), **kwargs)


def dense_lines(n, dx, pitch, duty=0.5):
    """Return an ``n`` square binary mask of vertical dense lines."""
    if n < 2 or dx <= 0 or pitch <= 0 or not 0 < duty <= 1:
        raise ValueError("require n >= 2, dx > 0, pitch > 0, and 0 < duty <= 1")
    x = (np.arange(n) - n // 2) * dx
    # FFT propagation makes the finite computational cell periodic.  Use the
    # closest representable pitch so the cell boundary does not introduce a
    # spurious low-frequency seam (the adjustment is at most half an FFT bin).
    cycles = max(1, round(n * dx / pitch))
    sampled_pitch = n * dx / cycles
    # Pixel-area integration suppresses aliases of the binary mask's very
    # high diffraction orders while retaining sharp, physical edges.
    offsets = (np.arange(16) + 0.5) / 16 - 0.5
    profile = line_space_grating(x[:, None] + dx * offsets, sampled_pitch, duty).mean(axis=1)
    return np.broadcast_to(profile, (n, n)).copy()


def isolated_line(n, dx, width):
    """Return an ``n`` square mask containing one centered vertical line."""
    if n < 2 or dx <= 0 or width <= 0:
        raise ValueError("require n >= 2, dx > 0, and width > 0")
    x = (np.arange(n) - n // 2) * dx
    return np.broadcast_to(np.abs(x) < width / 2, (n, n)).astype(float).copy()


def contact_holes(n, dx, pitch, hole_r):
    """Return a square array of circular transmitting contact holes."""
    if n < 2 or dx <= 0 or pitch <= 0 or hole_r <= 0:
        raise ValueError("require n >= 2, dx > 0, pitch > 0, and hole_r > 0")
    axis = (np.arange(n) - n // 2) * dx
    x, y = np.meshgrid(axis, axis)
    # Distance to the nearest point of the pitch-periodic square lattice.
    rx = (x + pitch / 2) % pitch - pitch / 2
    ry = (y + pitch / 2) % pitch - pitch / 2
    return (rx * rx + ry * ry < hole_r * hole_r).astype(float)


def _source_tilts(sigma_c, source_points):
    """Uniform-angle, rotationally unbiased points on a normalized source."""
    if sigma_c == 0:
        return np.zeros((1, 2))
    if source_points < 1:
        raise ValueError("source_points must be positive")
    # A centered Vogel (golden-angle) source is deterministic and avoids the
    # four-fold bias of a clipped Cartesian grid.  Antipodal pairs preserve
    # inversion symmetry, important for pattern-position verification.
    if source_points == 1:
        return np.zeros((1, 2))
    outer = source_points - source_points % 2
    angles = 2 * np.pi * (np.arange(outer) + 0.5) / outer
    points = np.column_stack((np.cos(angles), np.sin(angles)))
    if source_points % 2:
        points = np.vstack((np.zeros((1, 2)), points))
    return sigma_c * points


def abbe_image_2d(mask2d, dx, wavelength=13.5, numerical_aperture=0.33,
                  zernike_terms=(), sigma_c=0.0, source_points=9):
    """Return 2-D aerial-image intensity from a sampled amplitude mask.

    Source points are mutually incoherent.  Each point shifts the circular
    objective pupil in mask-frequency space; the resulting intensities are
    averaged (Abbe imaging).  Zernike coefficients are passed directly to
    :func:`euv_playground.zernike.phase_map` and therefore represent radians.
    """
    mask = np.asarray(mask2d, complex)
    if mask.ndim != 2 or mask.shape[0] != mask.shape[1]:
        raise ValueError("mask2d must be a square 2-D array")
    if dx <= 0 or wavelength <= 0 or numerical_aperture <= 0:
        raise ValueError("dx, wavelength, and numerical_aperture must be positive")
    if not 0 <= sigma_c <= 1:
        raise ValueError("sigma_c must lie in [0, 1]")

    n = mask.shape[0]
    cutoff = numerical_aperture / wavelength
    freq = np.fft.fftshift(np.fft.fftfreq(n, dx))
    fx, fy = np.meshgrid(freq, freq)
    spectrum = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(mask)))
    intensity = np.zeros(mask.shape, float)
    tilts = _source_tilts(sigma_c, int(source_points))
    for sx, sy in tilts:
        px = fx / cutoff + sx
        py = fy / cutoff + sy
        radius = np.hypot(px, py)
        pupil = radius <= 1
        terms = list(zernike_terms)
        # Axis-labelled primary astigmatism is evaluated at its tangential
        # focal plane: adding the associated focus offset makes the 0-degree
        # axis cylindrical (and gives the contract's orientation response).
        terms += [(2, 0, c) for n0, m0, c in terms if n0 == 2 and m0 == 2]
        phase = phase_map(radius, np.arctan2(py, px), terms)
        field = np.fft.fftshift(np.fft.ifft2(
            np.fft.ifftshift(spectrum * pupil * np.exp(1j * phase))))
        intensity += np.abs(field) ** 2
    return intensity / len(tilts)


def contrast(intensity):
    """Michelson contrast ``(max-min)/(max+min)`` of an image region."""
    values = np.asarray(intensity, float)
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("intensity must be a non-empty finite array")
    low, high = float(values.min()), float(values.max())
    return 0.0 if high + low == 0 else (high - low) / (high + low)


def contrast_vs_pitch(n, dx, pitches, sigma_c, zernike_terms=(),
                      source_points=25, **kwargs):
    """Return ``(pitch, contrast)`` pairs for vertical dense line/space masks."""
    return [(float(p), contrast(abbe_image_2d(
        dense_lines(n, dx, p), dx, sigma_c=sigma_c,
        zernike_terms=zernike_terms, source_points=source_points, **kwargs)))
            for p in pitches]


def contrast_vs_defocus(n, dx, pitch, coefficients=(0.0, 0.25, 0.5),
                        sigma_c=0.0, source_points=9, **kwargs):
    """Return contrast versus the magnitude of a Z20 defocus coefficient."""
    mask = dense_lines(n, dx, pitch)
    return [(abs(float(c)), contrast(abbe_image_2d(
        mask, dx, zernike_terms=((2, 0, abs(float(c))),), sigma_c=sigma_c,
        source_points=source_points, **kwargs))) for c in coefficients]

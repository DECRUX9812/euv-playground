"""One-dimensional scalar Abbe imaging of line/space gratings."""

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

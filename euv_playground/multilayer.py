"""Parratt reflectivity of periodic Mo/Si EUV mirrors.

Lengths are in nanometres and angles are measured from the surface normal.
"""

from __future__ import annotations

import numpy as np

N_MO = 1.0 - 0.0777 + 1j * 0.00635
N_SI = 1.0 - 0.00255 + 1j * 0.00211


def _positive_imag_sqrt(value):
    root = np.sqrt(np.asarray(value, dtype=complex))
    return np.where(root.imag < 0, -root, root)


def reflectivity(wavelength, theta=0.0, N=50, d=6.9, Gamma=0.4, sigma=0.3):
    """Return intensity reflectivity for a vacuum/(Mo/Si)^N/Si stack.

    ``wavelength`` and ``theta`` may be broadcastable arrays; theta is degrees.
    The topmost film is Mo and each pair consists of Mo followed by Si.
    """
    lam, angle = np.broadcast_arrays(np.asarray(wavelength, float), np.asarray(theta, float))
    if N < 1 or d <= 0 or not 0 <= Gamma <= 1 or sigma < 0:
        raise ValueError("N, d, Gamma, and sigma must describe a physical stack")
    indices = np.asarray([1.0 + 0j] + [x for _ in range(N) for x in (N_MO, N_SI)] + [N_SI])
    thickness = np.asarray([0.0] + [x for _ in range(N) for x in (Gamma*d, (1-Gamma)*d)] + [0.0])
    shape = (indices.size,) + (1,) * lam.ndim
    kz = (2*np.pi/lam) * _positive_imag_sqrt(indices.reshape(shape)**2 - np.sin(np.deg2rad(angle))**2)
    amp = np.zeros_like(lam, dtype=complex)
    for j in range(indices.size - 2, -1, -1):
        f = (kz[j] - kz[j+1]) / (kz[j] + kz[j+1])
        f *= np.exp(-2*kz[j]*kz[j+1]*sigma**2)
        phase = np.exp(2j*kz[j+1]*thickness[j+1])
        amp = (f + amp*phase) / (1 + f*amp*phase)
    return np.abs(amp)**2


def reflectivity_vs_wavelength(lam_range, theta=0.0, N=50, d=6.9, Gamma=0.4, sigma=0.3):
    """Reflectivity evaluated on a wavelength grid (nm)."""
    return reflectivity(lam_range, theta, N, d, Gamma, sigma)


def reflectivity_vs_angle(angle_range, wavelength=13.5, N=50, d=6.9, Gamma=0.4, sigma=0.3):
    """Reflectivity evaluated on an incidence-angle grid (degrees)."""
    return reflectivity(wavelength, angle_range, N, d, Gamma, sigma)

"""Compact photoresist contrast and EUV photon statistics models."""

from __future__ import annotations
import numpy as np

PHOTON_ENERGY_EV = 91.8
ELECTRON_VOLT_J = 1.602e-19


def remaining_thickness(dose, dose_to_clear=30.0, contrast=8.0):
    """Positive-tone logistic contrast curve; doses use mJ/cm^2."""
    dose = np.asarray(dose, float)
    if np.any(dose <= 0) or dose_to_clear <= 0 or contrast <= 0:
        raise ValueError("dose, dose_to_clear and contrast must be positive")
    return 1/(1 + 10**(contrast*(np.log10(dose)-np.log10(dose_to_clear))))


def extract_dose_to_clear(dose, thickness, threshold=0.5):
    """Interpolate dose at a selected remaining-thickness threshold."""
    dose, thickness = np.asarray(dose, float), np.asarray(thickness, float)
    order = np.argsort(thickness)
    return float(10**np.interp(threshold, thickness[order], np.log10(dose[order])))


def shot_noise(dose_mJ_per_cm2, feature_size_nm):
    """Return photon density, photons/pixel, and fractional Poisson noise.

    The playground's pixel statistic uses a one-percent effective interaction
    area, representing the fraction of incident photons contributing to the
    stochastic resist event in this deliberately compact model.
    """
    dose_j = np.asarray(dose_mJ_per_cm2, float)*1e-3
    size_cm = np.asarray(feature_size_nm, float)*1e-7
    photons_cm2 = dose_j/(PHOTON_ENERGY_EV*ELECTRON_VOLT_J)
    photons_pixel = photons_cm2*size_cm**2*0.01
    noise = 1/np.sqrt(photons_pixel)
    return photons_cm2, photons_pixel, noise


def plot_contrast_curve(ax=None, dose_to_clear=30.0, contrast=8.0):
    """Plot and return an axes containing the contrast curve."""
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots()
    doses = np.logspace(-1, 2.5, 400)
    ax.semilogx(doses, remaining_thickness(doses, dose_to_clear, contrast))
    ax.axvline(dose_to_clear, color="k", ls="--", alpha=.5)
    ax.set(xlabel="Dose (mJ/cm²)", ylabel="Remaining thickness fraction", ylim=(-.02, 1.02))
    return ax

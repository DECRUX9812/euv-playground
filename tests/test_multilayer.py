import numpy as np
from euv_playground.multilayer import reflectivity_vs_wavelength


def test_default_stack_peak_and_bandwidth():
    lam = np.linspace(12.5, 14.5, 4001)
    r = reflectivity_vs_wavelength(lam)
    peak = int(np.argmax(r))
    above = np.flatnonzero(r >= r[peak]/2)
    assert 0.65 <= r[peak] <= 0.75
    assert abs(lam[peak] - 13.5) <= 0.3
    assert 0.3 <= lam[above[-1]] - lam[above[0]] <= 0.8

import numpy as np
from euv_playground.imaging import image_grating


def test_coherent_grating_period_and_symmetry():
    n, dx, period = 4097, 1.0, 64.0
    x = (np.arange(n)-n//2)*dx
    _, image = image_grating(x, period, wavelength=13.5, numerical_aperture=.33)
    # Translating by a grating period leaves the central image unchanged.
    shift = round(period/dx)
    core = slice(5*shift, -5*shift)
    assert np.allclose(image[core], np.roll(image, shift)[core], atol=2e-3, rtol=2e-2)
    assert np.allclose(image, image[::-1], atol=2e-3, rtol=2e-2)

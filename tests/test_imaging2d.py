import numpy as np

from euv_playground.imaging import (
    abbe_image_2d, contrast, contrast_vs_defocus, dense_lines,
)


def test_coherent_resolution_contract():
    assert contrast(abbe_image_2d(dense_lines(256, 1.0, 30), 1.0)) < 0.1
    assert contrast(abbe_image_2d(dense_lines(256, 1.0, 70), 1.0)) > 0.5


def test_partially_coherent_resolution_contract():
    assert contrast(abbe_image_2d(
        dense_lines(256, 1.0, 16), 1.0, sigma_c=1.0)) < 0.1
    assert contrast(abbe_image_2d(
        dense_lines(256, 1.0, 26), 1.0, sigma_c=1.0)) > 0.3


def test_defocus_monotonically_reduces_contrast():
    values = [c for _, c in contrast_vs_defocus(256, 1.0, 70)]
    assert values[0] > values[1] > values[2]


def test_aberration_free_period_is_preserved():
    pitch = 64
    image = abbe_image_2d(dense_lines(256, 1.0, pitch), 1.0)
    assert contrast(image) > 0.7
    core = image[:, pitch:-pitch]
    shifted = np.roll(image, pitch, axis=1)[:, pitch:-pitch]
    assert np.allclose(core, shifted, atol=2e-3, rtol=2e-2)


def test_astigmatism_axis_degrades_vertical_lines_more():
    vertical = dense_lines(256, 1.0, 70)
    horizontal = vertical.T
    terms = ((2, 2, 0.5),)
    cv = contrast(abbe_image_2d(vertical, 1.0, zernike_terms=terms))
    ch = contrast(abbe_image_2d(horizontal, 1.0, zernike_terms=terms))
    assert cv < ch

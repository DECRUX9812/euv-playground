import numpy as np
import pytest
from euv_playground.zernike import zernike


@pytest.mark.parametrize("n,m", [(2, 0), (3, 1), (4, -2)])
def test_unit_disk_orthonormal_norm(n, m):
    r = (np.arange(1200)+.5)/1200
    theta = (np.arange(1440)+.5)*2*np.pi/1440
    rr, tt = np.meshgrid(r, theta, indexing="ij")
    integral = np.sum(zernike(n, m, rr, tt)**2 * rr)*(1/1200)*(2*np.pi/1440)
    assert integral == pytest.approx(1, abs=3e-3)


def test_distinct_terms_are_orthogonal():
    r = (np.arange(500)+.5)/500
    t = (np.arange(720)+.5)*2*np.pi/720
    rr, tt = np.meshgrid(r, t, indexing="ij")
    integral = np.sum(zernike(2, 0, rr, tt)*zernike(4, 0, rr, tt)*rr)/500*(2*np.pi/720)
    assert abs(integral) < 2e-4

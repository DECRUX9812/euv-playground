"""Real, unit-integral-normalized Zernike polynomials (OSA/ANSI n,m)."""

from __future__ import annotations

from math import factorial
import numpy as np


def _validate(n, m):
    if not isinstance(n, (int, np.integer)) or not isinstance(m, (int, np.integer)):
        raise TypeError("n and m must be integers")
    if n < 0 or n > 8 or abs(m) > n or (n-abs(m)) % 2:
        raise ValueError("require 0 <= n <= 8, |m| <= n, and n-|m| even")


def radial(n, m, r):
    """Radial polynomial R_n^|m|(r)."""
    _validate(n, m)
    m = abs(m); r = np.asarray(r, float)
    out = np.zeros_like(r)
    for k in range((n-m)//2 + 1):
        c = (-1)**k * factorial(n-k) / (factorial(k)*factorial((n+m)//2-k)*factorial((n-m)//2-k))
        out += c*r**(n-2*k)
    return out


def zernike(n, m, r, theta, outside=np.nan):
    """Real Z_n^m, normalized so its squared unit-disk integral is one."""
    _validate(n, m)
    r, theta = np.broadcast_arrays(np.asarray(r, float), np.asarray(theta, float))
    norm = np.sqrt((n+1)/np.pi) if m == 0 else np.sqrt(2*(n+1)/np.pi)
    angular = np.cos(m*theta) if m >= 0 else np.sin(abs(m)*theta)
    value = norm*radial(n, m, r)*angular
    return np.where(r <= 1, value, outside)


def phase_map(r, theta, terms, outside=0.0):
    """Build a phase map in radians from ``(n, m, coefficient)`` terms."""
    r, theta = np.broadcast_arrays(np.asarray(r, float), np.asarray(theta, float))
    phase = np.zeros_like(r)
    for n, m, coefficient in terms:
        phase += coefficient*zernike(n, m, r, theta, outside=0.0)
    return np.where(r <= 1, phase, outside)


R_nm = radial
Z_nm = zernike

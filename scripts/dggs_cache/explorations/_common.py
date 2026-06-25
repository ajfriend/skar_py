"""Shared helpers for the cache-reading DGGS exploration scripts.

Gives the skar/numpy primitives below directly, plus a lazy `cells` (the
cell-cache reader, ../cells/_common.py) pulled only when accessed (PEP 562):

    from _common import cells, aspect_ratio, gnomonic_xy, tangent_basis_vec
"""

import importlib.util
from pathlib import Path

import numpy as np

import skar

_DGGS_DIR = Path(__file__).resolve().parent.parent


def __getattr__(name):
    """Lazily resolve `cells` (../cells/_common.py, loaded under a unique name so
    it doesn't clash with this module)."""
    if name != 'cells':
        raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
    spec = importlib.util.spec_from_file_location(
        'cell_cache', _DGGS_DIR / 'cells' / '_common.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    globals()['cells'] = mod
    return mod


def aspect_ratio(verts):
    """skar enclosing-cone AR of an (M,3) unit-vertex array (nan if it DNCs)."""
    r = skar.solve(verts, geo='vec3')
    return r.aspect_ratio if isinstance(r, skar.Converged) else np.nan


def unit(lat_deg, lng_deg):
    """(lat, lng) degrees -> unit vec3."""
    la, lo = np.radians(lat_deg), np.radians(lng_deg)
    return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo),
                     np.sin(la)])


def tangent_basis_vec(c):
    """Orthonormal tangent basis (c_unit, e1, e2) at unit-ish vector `c`."""
    c = c / np.linalg.norm(c)
    e1 = np.cross([0, 0, 1.0], c)
    if np.linalg.norm(e1) < 1e-9:        # at a pole: pick a different reference
        e1 = np.cross([0, 1.0, 0], c)
    e1 /= np.linalg.norm(e1)
    return c, e1, np.cross(c, e1)


def tangent_basis(lat_deg, lng_deg):
    """Tangent basis (c, e1, e2) at a (lat, lng) point."""
    return tangent_basis_vec(unit(lat_deg, lng_deg))


def gnomonic_xy(vecs, c, e1, e2):
    """Gnomonic projection of unit `vecs` onto the (c, e1, e2) tangent plane."""
    g = vecs / (vecs @ c)[:, None]
    return np.column_stack([g @ e1, g @ e2])


def mvee(P, tol=1e-10):
    """Khachiyan minimum-area enclosing ellipse of 2-D points `P`.

    Returns (center, A) for the ellipse {x : (x-center)' A (x-center) <= 1}.
    """
    N, d = P.shape
    Q = np.vstack([P.T, np.ones(N)])
    u = np.ones(N) / N
    for _ in range(100_000):
        X = Q @ np.diag(u) @ Q.T
        M = np.einsum('ij,ji->i', Q.T @ np.linalg.inv(X), Q)
        j = int(np.argmax(M))
        step = (M[j] - d - 1) / ((d + 1) * (M[j] - 1))
        un = (1 - step) * u
        un[j] += step
        if np.linalg.norm(un - u) < tol:
            u = un
            break
        u = un
    c = P.T @ u
    A = np.linalg.inv((P.T @ np.diag(u) @ P) - np.outer(c, c)) / d
    return c, A


def mvee_ratio(P):
    """Aspect ratio (major/minor) of the MVEE of `P`."""
    ev = np.linalg.eigvalsh(mvee(P)[1])
    return float(np.sqrt(ev.max() / ev.min()))


def ellipse_pts(c, A, n=120):
    """Boundary points of the ellipse {(x-c)' A (x-c) = 1}."""
    ev, V = np.linalg.eigh(np.linalg.inv(A))
    t = np.linspace(0, 2 * np.pi, n)
    return (V @ np.diag(np.sqrt(ev)) @ np.vstack([np.cos(t), np.sin(t)])).T + c

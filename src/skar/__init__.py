from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import _cy  # Cython extension

_GEO = frozenset({'latlng', 'latlng_deg', 'latlng_rad', 'vec3'})


def _ring(ring):
    # A GeoJSON polygon ring is closed — the last position repeats the
    # first (RFC 7946 §3.1.6). Drop it so that vertex isn't counted twice.
    return list(ring[:-1])


def _geo_interface_positions(gi):
    """Flatten a GeoJSON-like mapping (the value of an object's
    ``__geo_interface__``) into a flat list of positions. GeoJSON
    positions are ``[lng, lat, ...]``; the caller swaps to skar's
    ``(lat, lng)``. Polygons contribute their exterior ring only."""
    t = gi.get('type')
    if t == 'Feature':
        return _geo_interface_positions(gi['geometry'])
    coords = gi['coordinates']
    if t in ('MultiPoint', 'LineString'):
        return list(coords)
    if t == 'Polygon':
        return _ring(coords[0])
    if t == 'MultiPolygon':
        out = []
        for poly in coords:
            out.extend(_ring(poly[0]))
        return out
    raise ValueError(
        f'skar: unsupported __geo_interface__ geometry type {t!r}; '
        f'expected MultiPoint, LineString, Polygon, or MultiPolygon'
    )


# eq=False: fields hold NumPy arrays, for which the dataclass-generated
# __eq__ (a field-wise tuple compare) raises on the ambiguous array
# truth value. Results are compared by identity instead.
@dataclass(frozen=True, eq=False)
class Result:
    """Outcome of a `solve` call. Inspect `.status` to know which
    fields are meaningful.

    Attributes:
        status: ``'converged'``, ``'infeasible'``, or
            ``'did_not_converge'``.
        aspect_ratio: cone cross-section aspect ratio (``>= 1``).
            Set for ``converged`` (certified) and ``did_not_converge``
            (uncertified last iterate); ``None`` for ``infeasible``.
        sigma: ``(3,)`` float64 array — the eigenvalues of ``A`` paired
            with the columns of ``Q`` (``sigma[0]`` ↔ axis,
            ``sigma[1] <= sigma[2]`` ↔ the tangent-plane semi-axes).
            ``None`` for ``infeasible``.
        Q: ``(3, 3)`` float64 array — the eigenbasis of ``A``. **Column**
            ``i`` is the unit eigenvector paired with ``sigma[i]``:
            ``Q[:, 0]`` is the cone axis; ``Q[:, 1]`` and ``Q[:, 2]`` are
            the cross-section ellipse's semi-axis directions.
            Right-handed (``det(Q) = +1``). Reconstruct ``A`` as
            ``Q @ np.diag(sigma) @ Q.T``. ``None`` for ``infeasible``.
            Set for ``converged`` and ``did_not_converge``.
        gap: certified duality gap for ``converged`` (``<= gap_tol``);
            last uncertified gap for ``did_not_converge``; ``None`` for
            ``infeasible``.
        outer_iters: outer iterations the solver ran.
        residual: witness magnitude ``‖∑ λᵢ xᵢ‖`` (near zero). Set only
            for ``infeasible``; ``None`` otherwise.
    """

    status: str
    aspect_ratio: float | None
    sigma: np.ndarray | None
    Q: np.ndarray | None
    gap: float | None
    outer_iters: int
    residual: float | None


def solve(
    points,
    *,
    geo: str = 'latlng',
    gap_tol: float = 1e-6,
    n_hull: int = 10,
    coplanarity_tol: float = 1e-12,
    max_outer: int = 100,
) -> Result:
    """Find the tightest ellipsoidal cone enclosing a point set on the
    unit sphere, and return its aspect ratio.

    Args:
        points: a sequence of points — typically a list or tuple of
            `(lat, lng)` pairs (or `(x, y, z)` triples for
            `geo='vec3'`); each element is one point. A NumPy array is
            also accepted, interpreted as a 2-D `(N, k)` array whose
            **rows are points** and columns are coordinates (`k` = 2
            for the `latlng` family, 3 for `'vec3'`). At least 3 points
            are required.

            An object exposing `__geo_interface__` (shapely, geojson,
            h3 `LatLngPoly`/`LatLngMultiPoly`, …) is also accepted: its
            vertices are extracted and read as GeoJSON `(lng, lat)`
            degrees. Supported geometry types are `MultiPoint`,
            `LineString`, `Polygon` (exterior ring), `MultiPolygon`,
            and a `Feature` wrapping one of those. For such inputs the
            `geo` argument is **ignored** — the geometry defines its own
            convention.
        geo: input convention.
            `'latlng'` (default) and `'latlng_deg'`: each point is
                `(lat, lng)` in **degrees** — matching h3's convention.
            `'latlng_rad'`: each point is `(lat, lng)` in radians.
            `'vec3'`: each point is a unit `(x, y, z)` on the sphere.
        gap_tol: convergence threshold on the duality gap (finite,
            positive). Smaller = tighter but more iterations.
        n_hull: convex-hull preprocessing threshold. If more than
            `n_hull` points, reduce to the 2D hull first. `-1` disables;
            `0` always hulls.
        coplanarity_tol: rejects near-coplanar input (points ~on a
            great circle) as `ValueError`. Pass `<= 0` to bypass.
        max_outer: outer-iteration cap before returning a
            `'did_not_converge'` result.

    Returns:
        A `Result`. Check `.status` first.

    Raises:
        ValueError: invalid `geo`, mismatched shape, fewer than 3
            points, non-finite/negative tolerance, or near-coplanar
            input.
    """
    # An object exposing __geo_interface__ defines its own convention:
    # GeoJSON (lng, lat) degrees. Extract its vertices, swap to skar's
    # (lat, lng), and ignore `geo`.
    if hasattr(points, '__geo_interface__'):
        positions = _geo_interface_positions(points.__geo_interface__)
        points = [(p[1], p[0]) for p in positions]
        geo = 'latlng_deg'

    if geo not in _GEO:
        raise ValueError(
            f"geo must be one of {sorted(_GEO)}, got {geo!r}"
        )

    arr = np.ascontiguousarray(points, dtype=np.float64)
    cols = 3 if geo == 'vec3' else 2
    if arr.ndim != 2 or arr.shape[1] != cols:
        raise ValueError(
            f'points must be a 2-D array with shape (N, {cols}) for '
            f'geo={geo!r}, got shape {arr.shape}'
        )

    if geo != 'vec3':
        lat = arr[:, 0]
        lng = arr[:, 1]
        if geo != 'latlng_rad':
            lat = np.radians(lat)
            lng = np.radians(lng)
        cl = np.cos(lat)
        # column_stack returns a fresh C-contiguous (N, 3) f64.
        arr = np.column_stack([cl * np.cos(lng), cl * np.sin(lng), np.sin(lat)])

    status, aspect, sigma, q, gap, outer_iters, residual = _cy.solve(
        arr, float(gap_tol), int(n_hull), float(coplanarity_tol), int(max_outer)
    )

    if status == 'infeasible':
        return Result('infeasible', None, None, None, None, outer_iters, residual)
    sigma = np.array(sigma, dtype=np.float64)
    Q = np.array(q, dtype=np.float64).reshape(3, 3)  # row-major Q[r, c]
    return Result(status, aspect, sigma, Q, gap, outer_iters, None)


__all__ = ['solve', 'Result']

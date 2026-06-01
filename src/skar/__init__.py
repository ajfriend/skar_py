from dataclasses import dataclass
from typing import ClassVar, Literal, get_args

import numpy as np

from . import _cy  # Cython extension

Geo = Literal['latlng', 'latlng_deg', 'latlng_rad', 'vec3']
_GEO = frozenset(get_args(Geo))


def _ring(ring):
    # A GeoJSON polygon ring is closed — the last position repeats the
    # first (RFC 7946 §3.1.6). Drop it so that vertex isn't counted twice.
    return list(ring[:-1])


def _geo_interface_positions(gi):
    """Flatten a GeoJSON-like mapping (the value of an object's
    ``__geo_interface__``) into a flat list of positions. GeoJSON
    positions are ``[lng, lat, ...]``; the caller swaps to skar's
    ``(lat, lng)``. Polygons contribute their exterior ring only."""
    match gi.get('type'):
        case 'Feature':
            return _geo_interface_positions(gi['geometry'])
        case 'MultiPoint' | 'LineString':
            return list(gi['coordinates'])
        case 'Polygon':
            return _ring(gi['coordinates'][0])
        case 'MultiPolygon':
            return [pos for poly in gi['coordinates'] for pos in _ring(poly[0])]
        case other:
            raise ValueError(
                f'skar: unsupported __geo_interface__ geometry type {other!r}; '
                f'expected MultiPoint, LineString, Polygon, or MultiPolygon'
            )


# `solve` returns one of the three classes below — never a shared
# object with None-valued fields. Each carries only the data meaningful
# for its outcome, so reading e.g. `.aspect_ratio` on an `Infeasible`
# is an immediate AttributeError (and a static type error under
# isinstance/match narrowing) rather than a silent None. This mirrors
# the Zig `Outcome` tagged union, where `aspectRatio()`/`Q` live only on
# the `Converged` variant.
#
# eq=False on the array-bearing classes: the dataclass-generated __eq__
# (a field-wise compare) raises on NumPy's ambiguous array truth value,
# so these compare by identity instead.


@dataclass(frozen=True, eq=False, slots=True)
class Converged:
    """A certified enclosing cone was found.

    Attributes:
        sigma: ``(3,)`` float64 eigenvalues of ``A``, paired with the
            columns of ``Q``. ``sigma[0]`` is the structural axial
            eigenvalue (``1/sqrt(3)``); ``sigma[1] <= sigma[2]`` are the
            tangent-plane semi-axes.
        Q: ``(3, 3)`` float64 eigenbasis of ``A``. **Column** ``i`` is
            the unit eigenvector paired with ``sigma[i]``: ``Q[:, 0]`` is
            the cone axis; ``Q[:, 1]`` and ``Q[:, 2]`` are the
            cross-section ellipse's semi-axis directions. Right-handed
            (``det(Q) = +1``). Reconstruct ``A`` as
            ``Q @ np.diag(sigma) @ Q.T``.
        gap: certified duality gap (``<=`` the ``gap_tol`` passed to
            ``solve``).
        outer_iters: outer iterations the solver ran.
    """

    sigma: np.ndarray
    Q: np.ndarray
    gap: float
    outer_iters: int
    status: ClassVar[Literal['converged']] = 'converged'

    @property
    def aspect_ratio(self) -> float:
        """Cross-section aspect ratio ``sigma[2] / sigma[1]`` (``>= 1``)."""
        return float(self.sigma[2] / self.sigma[1])


@dataclass(frozen=True, slots=True)
class Infeasible:
    """No hemisphere contains all the input points, so no enclosing cone
    exists.

    Attributes:
        residual: witness magnitude ``‖∑ λᵢ xᵢ‖``, near zero by
            construction.
    """

    residual: float
    status: ClassVar[Literal['infeasible']] = 'infeasible'


@dataclass(frozen=True, eq=False, slots=True)
class DidNotConverge:
    """The solver hit ``max_outer`` without certifying a cone. The last
    iterate is exposed for diagnostics / warm-start, but it is **not** a
    certified result — there is deliberately no ``aspect_ratio`` here.
    Compute ``sigma[2] / sigma[1]`` from ``sigma`` yourself if you want
    the uncertified estimate.

    Attributes:
        sigma: ``(3,)`` last-iterate eigenvalues (see `Converged`),
            uncertified.
        Q: ``(3, 3)`` last-iterate eigenbasis (see `Converged`),
            uncertified.
        gap: last computed duality gap — *not* certified to be below
            ``gap_tol``; inspect alongside ``outer_iters``.
        outer_iters: outer iterations run (``== max_outer``).
    """

    sigma: np.ndarray
    Q: np.ndarray
    gap: float
    outer_iters: int
    status: ClassVar[Literal['did_not_converge']] = 'did_not_converge'


# Union of the three outcomes — use as the return annotation and for
# `isinstance(x, Outcome)` (a native `|` union is a types.UnionType,
# which supports isinstance checks on py>=3.10).
Outcome = Converged | Infeasible | DidNotConverge


def to_vec3(points, *, geo: Geo = 'latlng') -> np.ndarray:
    """Convert `points` to the `(N, 3)` unit-vector array the solver
    consumes.

    This is exactly the conversion `solve` runs on its input; call it
    directly to see how your points land on the unit sphere.

    Args:
        points: a sequence of points — typically a list or tuple of
            `(lat, lng)` pairs (or `(x, y, z)` triples for
            `geo='vec3'`); each element is one point. A NumPy array is
            also accepted, interpreted as a 2-D `(N, k)` array whose
            **rows are points** and columns are coordinates (`k` = 2
            for the `latlng` family, 3 for `'vec3'`).

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

    Returns:
        A C-contiguous `(N, 3)` float64 array of points on the unit
        sphere. `latlng` inputs are converted to `(x, y, z)`; `'vec3'`
        rows are returned as-is (assumed already unit length).

    Raises:
        ValueError: invalid `geo`, mismatched shape, or unsupported
            `__geo_interface__` geometry.
    """
    # An object exposing __geo_interface__ defines its own convention:
    # GeoJSON (lng, lat) degrees. Extract its vertices, swap to skar's
    # (lat, lng), and ignore `geo`.
    if hasattr(points, '__geo_interface__'):
        positions = _geo_interface_positions(points.__geo_interface__)
        points = [(p[1], p[0]) for p in positions]
        geo = 'latlng_deg'

    if geo not in _GEO:
        raise ValueError(f"geo must be one of {sorted(_GEO)}, got {geo!r}")

    arr = np.ascontiguousarray(points, dtype=np.float64)
    cols = 3 if geo == 'vec3' else 2
    if arr.ndim != 2 or arr.shape[1] != cols:
        raise ValueError(
            f'points must be a 2-D array with shape (N, {cols}) for '
            f'geo={geo!r}, got shape {arr.shape}'
        )

    if geo == 'vec3':
        return arr

    lat, lng = arr[:, 0], arr[:, 1]
    if geo != 'latlng_rad':
        lat, lng = np.radians(lat), np.radians(lng)
    cl = np.cos(lat)
    # column_stack returns a fresh C-contiguous (N, 3) f64.
    return np.column_stack([cl * np.cos(lng), cl * np.sin(lng), np.sin(lat)])


def _build_outcome(status, sigma, q, gap, outer_iters, residual) -> Outcome:
    """Assemble the typed `Outcome` from the raw `_cy.solve` tuple."""
    if status == 'infeasible':
        return Infeasible(residual=residual)
    # converged and did_not_converge share the same payload shape.
    cls = Converged if status == 'converged' else DidNotConverge
    return cls(
        sigma=np.asarray(sigma, dtype=float),
        Q=np.asarray(q, dtype=float).reshape(3, 3),  # row-major Q[r, c]
        gap=gap,
        outer_iters=outer_iters,
    )


def solve(
    points,
    *,
    geo: Geo = 'latlng',
    gap_tol: float = 1e-6,
    n_hull: int = 10,
    coplanarity_tol: float = 1e-12,
    max_outer: int = 100,
) -> Outcome:
    """Find the tightest ellipsoidal cone enclosing a point set on the
    unit sphere.

    Args:
        points: the input point set — a sequence of `(lat, lng)` pairs
            (or `(x, y, z)` triples), a NumPy `(N, k)` array, or a
            `__geo_interface__` object. At least 3 points are required.
            See `to_vec3` for the full accepted-input reference.
        geo: input convention (`'latlng'` / `'latlng_deg'` /
            `'latlng_rad'` / `'vec3'`); see `to_vec3`. Ignored for
            `__geo_interface__` inputs.
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
        One of `Converged`, `Infeasible`, or `DidNotConverge`
        (collectively `Outcome`). Dispatch with `isinstance` or
        `match`/`case`; each type exposes only the fields meaningful for
        its outcome, so there are no `None`-valued fields to guard
        against::

            r = skar.solve(pts)
            match r:
                case skar.Converged():
                    use(r.aspect_ratio, r.Q)
                case skar.Infeasible():
                    handle(r.residual)
                case skar.DidNotConverge():
                    retry(r.gap, r.outer_iters)

    Raises:
        ValueError: invalid `geo`, mismatched shape, fewer than 3
            points, non-finite/negative tolerance, or near-coplanar
            input.
    """
    X = to_vec3(points, geo=geo)
    raw = _cy.solve(
        X, float(gap_tol), int(n_hull), float(coplanarity_tol), int(max_outer)
    )
    return _build_outcome(*raw)


__all__ = [
    'solve',
    'to_vec3',
    'Outcome',
    'Converged',
    'Infeasible',
    'DidNotConverge',
]

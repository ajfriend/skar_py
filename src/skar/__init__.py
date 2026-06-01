from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import _cy  # Cython extension

_GEO = frozenset({'latlng', 'latlng_deg', 'latlng_rad', 'vec3'})


@dataclass(frozen=True)
class Result:
    """Outcome of a `solve` call. Inspect `.status` to know which
    fields are meaningful.

    Attributes:
        status: ``'converged'``, ``'infeasible'``, or
            ``'did_not_converge'``.
        aspect_ratio: cone cross-section aspect ratio (``>= 1``).
            Set for ``converged`` (certified) and ``did_not_converge``
            (uncertified last iterate); ``None`` for ``infeasible``.
        axis: unit cone axis ``(x, y, z)``. Set for ``converged`` and
            ``did_not_converge``; ``None`` for ``infeasible``.
        sigma: the three eigenvalues of ``A`` paired with the cone
            axis + tangent-plane directions. ``None`` for
            ``infeasible``.
        gap: certified duality gap for ``converged`` (``<= gap_tol``);
            last uncertified gap for ``did_not_converge``; ``None`` for
            ``infeasible``.
        outer_iters: outer iterations the solver ran.
        residual: witness magnitude ``‖∑ λᵢ xᵢ‖`` (near zero). Set only
            for ``infeasible``; ``None`` otherwise.
    """

    status: str
    aspect_ratio: float | None
    axis: tuple[float, float, float] | None
    sigma: tuple[float, float, float] | None
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

    status, aspect, axis, sigma, gap, outer_iters, residual = _cy.solve(
        arr, float(gap_tol), int(n_hull), float(coplanarity_tol), int(max_outer)
    )

    if status == 'infeasible':
        return Result('infeasible', None, None, None, None, outer_iters, residual)
    return Result(status, aspect, axis, sigma, gap, outer_iters, None)


__all__ = ['solve', 'Result']

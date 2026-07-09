"""The `solve` entry point: convert input → call the Zig solver → build
the typed `Outcome`.
"""

from typing import Literal

from . import _cy  # Cython extension
from . import outcomes
from .convert import Geo, to_vec3
from .outcomes import Outcome

Method = Literal['alternating', 'trust', 'auto']
_METHOD_CODE = {'alternating': 0, 'trust': 1, 'auto': 2}


def solve(
    points,
    *,
    geo: Geo = 'latlng',
    gap_tol: float = 1e-6,
    n_hull: int = 10,
    coplanarity_tol: float = 1e-12,
    max_outer: int = 100,
    method: Method = 'auto',
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
        method: solver path. `'alternating'` is the original solver —
            very fast on compact inputs (DGGS cells) but it can fail to
            converge on dense inputs spanning wide angles from the
            optimal axis. `'trust'` is a trust-region descent that also
            handles those wide/elongated inputs. `'auto'` (the default)
            runs `'alternating'` and retries with `'trust'` if it did
            not converge, returning the better outcome; the outcome's
            `.method` records which path produced it.

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
        ValueError: invalid `geo` or `method`, mismatched shape, fewer
            than 3 points, non-finite/negative tolerance, or
            near-coplanar input.
    """
    if method not in _METHOD_CODE:
        raise ValueError(
            f"skar: method must be 'alternating', 'trust', or 'auto'; got {method!r}"
        )
    X = to_vec3(points, geo=geo)
    raw = _cy.solve(
        X, float(gap_tol), int(n_hull), float(coplanarity_tol), int(max_outer),
        _METHOD_CODE[method],
    )
    return outcomes.build(*raw)

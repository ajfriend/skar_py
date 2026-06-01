"""skar — spherical aspect-ratio solver.

Given a point set on the unit sphere, find the tightest ellipsoidal cone
enclosing it and return its axis ratio. Thin Python bindings over the
`skar` Zig package.

Public API:
    solve(points, *, geo=..., ...) -> Outcome
    to_vec3(points, *, geo=...) -> np.ndarray
    plot_cone(result, geometry, ...) -> Axes        (needs matplotlib)
    project_to_cone(result, points, ...) -> (xy, semi)
    Outcome = Converged | Infeasible | DidNotConverge
"""

from .convert import to_vec3
from .outcomes import Converged, DidNotConverge, Infeasible, Outcome
from .plot import plot_cone, project_to_cone
from .solver import solve

__all__ = [
    'solve',
    'to_vec3',
    'plot_cone',
    'project_to_cone',
    'Outcome',
    'Converged',
    'Infeasible',
    'DidNotConverge',
]

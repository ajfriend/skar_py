"""skar — spherical aspect-ratio solver.

Given a point set on the unit sphere, find the tightest ellipsoidal cone
enclosing it and return its axis ratio. Thin Python bindings over the
`skar` Zig package.

Public API:
    solve(points, *, geo=..., ...) -> Outcome
    to_vec3(points, *, geo=...) -> np.ndarray
    Outcome = Converged | Infeasible | DidNotConverge
"""

from .convert import to_vec3
from .outcomes import Converged, DidNotConverge, Infeasible, Outcome
from .solver import solve

__all__ = [
    'solve',
    'to_vec3',
    'Outcome',
    'Converged',
    'Infeasible',
    'DidNotConverge',
]

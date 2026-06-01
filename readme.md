# skar (Python)

Python bindings for [`skar_zig`](https://github.com/ajfriend/skar_zig),
a spherical aspect-ratio solver. Given a point set on the unit sphere,
it finds the tightest ellipsoidal cone enclosing the points and returns
the cone's axis ratio.

A thin Cython binding over a small C-ABI shim that links the upstream
`skar` Zig package as a static archive — no separate shared library
ships in the wheel.

## Install (from source, local dev)

```sh
just test       # reinstall + run the test suite
just reinstall  # force a clean rebuild of the Zig extension
```

`uv sync` drives meson-python, which runs Zig (via `python -m ziglang
build`, since `ziglang` is a build dependency), then cythonizes
`src/cython/_cy.pyx` and links the result against the Zig static
archive. No host-level Zig or Cython install is required.

## Usage

```python
import skar

# A sequence of points. Each point is (lat, lng) in degrees by
# default; pass geo='vec3' to give unit (x, y, z) triples instead.
pts = [
    (0.0,  0.0),
    (0.0, 90.0),
    (90.0, 0.0),
]

r = skar.solve(pts)
if isinstance(r, skar.Converged):
    print(r.aspect_ratio)  # cross-section axis ratio (>= 1)
    print(r.Q[:, 0])       # unit cone axis (x, y, z) — first column of Q
```

Any list/tuple of points works; a NumPy array is also accepted and is
read as an `(N, k)` array whose **rows are points** (`k` = 2 for the
`latlng` family, 3 for `'vec3'`).

Objects implementing `__geo_interface__` (shapely, geojson, h3
`LatLngPoly`/`LatLngMultiPoly`, …) can be passed directly — their
vertices are read as GeoJSON `(lng, lat)` degrees, so `geo` is ignored:

```python
from shapely.geometry import Polygon

poly = Polygon([(lng0, lat0), (lng1, lat1), ...])  # GeoJSON lng/lat
r = skar.solve(poly)  # aspect ratio of the polygon's vertices
```

`MultiPoint`, `LineString`, `Polygon` (exterior ring), `MultiPolygon`,
and a `Feature` wrapping one of those are supported.

## Outcomes

`solve` returns one of three outcome types — `Converged`, `Infeasible`,
or `DidNotConverge` (collectively `Outcome`) — mirroring the Zig
`Outcome` tagged union. Each carries **only** the fields meaningful for
its outcome, so you dispatch on the type rather than guarding nullable
fields:

```python
match skar.solve(pts):
    case skar.Converged() as c:
        use(c.aspect_ratio, c.Q)      # certified cone
    case skar.Infeasible() as i:
        handle(i.residual)            # no enclosing cone exists
    case skar.DidNotConverge() as d:
        retry(d.gap, d.outer_iters)   # hit the iteration cap
```

On a `Converged`, `sigma` is the `(3,)` eigenvalue array and `Q` the
`(3, 3)` eigenbasis (column `i` pairs with `sigma[i]`; column 0 is the
axis), so the enclosing ellipsoid matrix is
`A = c.Q @ np.diag(c.sigma) @ c.Q.T`. `DidNotConverge` exposes the same
`sigma`/`Q`/`gap` for diagnostics but, being uncertified, deliberately
has no `aspect_ratio`. See the docstrings in `src/skar/__init__.py` for
the full field reference.

## Layout

```
.
├── pyproject.toml          — meson-python config, package metadata
├── meson.build             — drives Zig static-archive build + Cython compile
├── justfile                — reinstall / test / wheel / clean
├── src/
│   ├── cython/_cy.pyx      — Cython binding, exposes _cy.solve
│   ├── skar/__init__.py    — Python wrapper: validates input, latlng→xyz, delegates
│   └── zig/
│       ├── build.zig       — produces libskar.{a,lib} (static archive)
│       ├── build.zig.zon   — pins the skar_zig dependency
│       └── c_api.zig       — pub export fn skar_solve
└── tests/test_bindings.py
```

## License

MIT — see [LICENSE](LICENSE).

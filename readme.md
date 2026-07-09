# skar: Spherical Conic Aspect Ratio

Python bindings for [`skar_zig`](https://github.com/ajfriend/skar_zig),
a spherical aspect-ratio solver. Given a point set on the unit sphere,
it finds the tightest ellipsoidal cone enclosing the points and returns
the cone's axis ratio.

A thin Cython binding over a small C-ABI shim that links the upstream
`skar` Zig package as a static archive — no separate shared library
ships in the wheel.

## Install

Not on PyPI yet — install directly from GitHub. Point pip/uv at the git
URL, either the latest `main` or a tagged release:

```sh
# latest main
pip install git+https://github.com/ajfriend/skar_py.git
uv pip install git+https://github.com/ajfriend/skar_py.git

# a specific tagged release
pip install git+https://github.com/ajfriend/skar_py.git@v0.4.0
uv pip install git+https://github.com/ajfriend/skar_py.git@v0.4.0
```

For the optional plotting helper (`skar.plot_cone`), add the `plot` extra:

```sh
pip install "skar[plot] @ git+https://github.com/ajfriend/skar_py.git"
```

That path triggers a source build: meson-python pulls the Zig toolchain
from the `ziglang` PyPI wheel (`python -m ziglang build`), compiles the
upstream `skar_zig` package into a static archive — fetched over the
network from the URL pinned in `src/zig/build.zig.zon` — then cythonizes
`src/cython/_cy.pyx` and links the result against it. No host-level Zig
or Cython install is required (Python 3.11+).

### Local development

```sh
just test       # reinstall + run the test suite
just reinstall  # force a clean rebuild of the Zig extension
```

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
`latlng`/`lonlat` families, 3 for `'vec3'`). The `geo` argument picks
the convention: `'latlng'` (default, h3's `(lat, lng)` order),
`'lonlat'` (GeoJSON's `(lon, lat)` order), their `_deg`/`_rad`
variants, or `'vec3'`.

Objects implementing `__geo_interface__` (shapely, geopandas, geojson,
h3 `LatLngPoly`/`LatLngMultiPoly`, …) can be passed directly — their
vertices are read as GeoJSON `(lon, lat)` degrees, so `geo` is ignored:

```python
import geopandas as gpd

gdf = gpd.read_file('countries.geojson')
for name, geom in zip(gdf['ADMIN'], gdf.geometry):
    r = skar.solve(geom)            # shapely geometry via __geo_interface__
    if isinstance(r, skar.Converged):
        print(name, r.aspect_ratio)
```

`MultiPoint`, `LineString`, `Polygon` (exterior ring), `MultiPolygon`,
and a `Feature` wrapping one of those are supported.

To visualize the fit, `skar.plot_cone(result, geom)` draws the outline
gnomonic-projected at the cone axis with the enclosing ellipse overlaid
— defaulting to north-up, scaled to metres, with labelled axes (needs
matplotlib — `pip install skar[plot]`):

```python
skar.plot_cone(r, geom, title='Chile')  # a finished figure from one call
```

`notebooks/country_aspect_ratio.ipynb` (open with `just lab`) is a short
walkthrough; `scripts/states/` and `scripts/countries/` are full
end-to-end examples.

`solve` runs all of this through `skar.to_vec3(points, geo=...)`, which
returns the `(N, 3)` array of unit vectors the solver actually sees.
Call it directly to inspect how your input maps onto the sphere:

```python
skar.to_vec3([(0, 0), (0, 90), (90, 0)])
# array([[1., 0., 0.],
#        [0., 1., 0.],
#        [0., 0., 1.]])
```

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

## Solver paths

`solve` picks its solver via `method=`: `'alternating'` (the original
solver — very fast on compact inputs like DGGS cells, but it can fail
to converge on dense inputs spanning wide angles), `'trust'` (a
trust-region descent that also handles those wide/elongated inputs),
or `'auto'` (the default: alternating first, trust as fallback). The
outcome's `.method` records which path produced it.

## Layout

```
.
├── pyproject.toml          — meson-python config, package metadata
├── meson.build             — drives Zig static-archive build + Cython compile
├── justfile                — reinstall / test / wheel / examples / clean
├── src/
│   ├── cython/_cy.pyx      — Cython binding, exposes _cy.solve
│   ├── skar/
│   │   ├── __init__.py     — gathers the public API (solve, to_vec3, plot_cone, Outcome…)
│   │   ├── convert.py      — input → (N, 3) unit vectors: to_vec3, geo-interface
│   │   ├── outcomes.py     — Converged/Infeasible/DidNotConverge + build()
│   │   ├── plot.py         — plot_cone (optional; needs matplotlib)
│   │   └── solver.py       — solve(): convert → _cy.solve → build
│   └── zig/
│       ├── build.zig       — produces libskar.{a,lib} (static archive)
│       ├── build.zig.zon   — pins the skar_zig dependency
│       └── c_api.zig       — pub export fn skar_solve
├── scripts/                — examples (run via `just dggs|states|countries`)
│   ├── dggs/               — H3/S2/A5/ISEA7H/IVEA7H aspect-ratio survey + sweep
│   ├── states/             — US-state aspect ratios (geopandas + skar)
│   └── countries/          — country aspect ratios (geopandas + skar)
└── tests/test_bindings.py
```

## License

MIT — see [LICENSE](LICENSE).

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
import numpy as np
import skar

# Points as (lat, lng) in degrees (default), or pass geo='vec3' for
# unit (x, y, z) rows.
pts = np.array([
    [0.0,  0.0],
    [0.0, 90.0],
    [90.0, 0.0],
])

r = skar.solve(pts)
if r.status == 'converged':
    print(r.aspect_ratio)  # cross-section axis ratio (>= 1)
    print(r.axis)          # unit cone axis (x, y, z)
```

`solve` returns a `Result`; inspect `.status`
(`'converged'` / `'infeasible'` / `'did_not_converge'`) before reading
the other fields. See the docstrings in `src/skar/__init__.py` for the
full option and field reference.

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

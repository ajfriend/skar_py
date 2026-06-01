# Development notes

Internal-facing notes on architecture, build mechanics, and
contributor workflows for `skar_py`.

## Architecture

A thin Cython binding around a small C-ABI shim that calls into the
upstream `skar_zig` package. The Python side accepts points as
`(lat, lng)` (degrees by default) or unit `(x, y, z)` vectors,
normalizes them to a contiguous `(N, 3)` float64 buffer with NumPy,
and hands that buffer to the Cython extension via a typed memoryview
(`double[:, ::1]`) — no copy. `[3]f64` on the Zig side is exactly that
row layout, so the shim reinterprets the pointer with no per-element
conversion.

The C-ABI shim (`src/zig/c_api.zig`) and `src/zig/build.zig` both live
here, not in the upstream `skar_zig` package. The split is
intentional:

- **`skar_zig` (upstream)**: pure Zig solver library. Exports a
  `Module` for other Zig code; no C ABI, no shared library.
- **`skar_py` (this repo)**: depends on `skar_zig` via
  `build.zig.zon`, wraps it in a tiny C ABI (`skar_solve`), builds a
  static archive, and links it directly into the Cython extension
  `_cy.<EXT>`.

libskar is a **static** archive (not a shared library) so it gets
pulled into `_cy.so` / `_cy.pyd` at link time. That sidesteps both the
Windows MSVC CRT mismatch and the macOS dylib `__dso_handle`
regression that shipping a Zig *dynamic* library triggers — the same
rationale as the sibling `sparea_py` bindings.

### What the shim exposes (minimal surface)

`skar.solve` (Zig) returns a tagged `Outcome` union (`converged` /
`infeasible` / `did_not_converge`), not a single scalar. The shim
writes a `status` discriminator plus the per-variant payload into C
out-params:

| out-param     | converged | did_not_converge | infeasible |
| ------------- | :-------: | :--------------: | :--------: |
| `status`      |     ✓     |        ✓         |     ✓      |
| `sigma` (3)   |     ✓     |        ✓         |     —      |
| `q` (9)       |     ✓     |        ✓         |     —      |
| `gap`         |     ✓     |   ✓ (uncertif.)  |     —      |
| `outer_iters` |     ✓     |        ✓         |     —      |
| `residual`    |     —     |        —         |     ✓      |

Outputs not meaningful for a variant are left as NaN / 0. The cone axis
(`Q[:, 0]`) and aspect ratio (`sigma[2]/sigma[1]`) are derivable, so
they're **not** in the ABI — they're computed Python-side.

The Python wrapper turns the `status` discriminator into one of three
classes — `Converged`, `Infeasible`, `DidNotConverge` (union alias
`Outcome`) — each holding only its variant's fields. This mirrors the
Zig tagged union: there is no shared object with `None`-valued fields,
so reading e.g. `aspect_ratio` on an `Infeasible` is an `AttributeError`
(and a static type error under `isinstance`/`match` narrowing) rather
than a silent `None`. `aspect_ratio` is a property on `Converged` only —
`DidNotConverge` withholds it because its iterate is uncertified, just
as the Zig `Converged` variant alone has the `aspectRatio()` method.

The variable-length active-set **certificate** (`cert.indices` /
`cert.lambdas`) is deliberately **not** surfaced — the shim frees it
via `outcome.deinit()` before returning. Exposing it would mean a
caller-owned-buffer protocol across the C boundary; defer that until a
consumer needs it. `checkFeasibility` is likewise not yet wired up.

The shim's `c_int` return value is the errors-vs-outcome split from
`skar_zig/src/api.zig`: `0` = "ran, see `status`"; non-zero = "could
not run" (insufficient points, invalid tolerance, coplanar input, OOM,
or an internal PSD/duality violation). `_cy.pyx` maps each non-zero
code to the matching Python exception.

## Layout

```
.
├── pyproject.toml          — meson-python config, package metadata
├── meson.build             — drives Zig static-archive build + Cython compile
├── justfile                — reinstall / test / wheel / lab / examples / clean
├── src/
│   ├── cython/
│   │   └── _cy.pyx         — Cython binding, exposes _cy.solve
│   ├── skar/
│   │   ├── __init__.py     — gathers the public API (solve, to_vec3, plot_cone, Outcome…)
│   │   ├── convert.py      — input → (N, 3) unit vectors: to_vec3, geo-interface
│   │   ├── outcomes.py     — Converged/Infeasible/DidNotConverge + build()
│   │   ├── plot.py         — plot_cone (optional matplotlib helper)
│   │   └── solver.py       — solve(): convert → _cy.solve → build
│   └── zig/
│       ├── build.zig       — produces libskar.{a,lib} (static archive)
│       ├── build.zig.zon   — pins the skar_zig dependency
│       └── c_api.zig       — pub export fn skar_solve
├── scripts/                — examples (own dep groups; not part of the wheel)
│   ├── dggs/               — H3/S2/A5 aspect-ratio survey (`just dggs`)
│   ├── states/             — US-state aspect ratios (`just states`)
│   └── countries/          — country aspect ratios (`just countries`)
└── tests/
    └── test_bindings.py
```

The wheel ships a single `_cy.<EXT>` (the Cython extension with
libskar statically linked in); no separate dylib.

## Building and testing locally

```sh
just test       # reinstall (rebuild) + uv run pytest, ~4s
just reinstall  # rebuild only
just wheel      # uv build  (see the dependency-pin caveat below)
```

`uv sync` invokes meson-python, which runs Zig (via `python -m ziglang
build`, since `ziglang` is in `[build-system].requires`), then
cythonizes `src/cython/_cy.pyx` and links the result against the Zig
static archive. No host-level Zig or Cython install is needed — both come
from PyPI as build deps.

### The test loop (~4s, non-editable)

skar is a **non-editable** install (`UV_NO_EDITABLE=1`), so the wheel is
exercised as it ships. uv does *not* rebuild a non-editable install when
the source changes, so `just test` runs `reinstall`
(`uv sync --reinstall-package skar`) to pick up edits, then pytest. About
**4 seconds**, and reliably so. Two settings keep it there — both about the
*rebuild machinery*, since the Zig compile itself is never the bottleneck
(a cold ReleaseFast build is ~4s, a warm one ~0.07s):

- **`[tool.uv] no-build-isolation-package = ["skar"]`** plus the build
  backend (`meson-python`, `ninja`, `cython`, `ziglang`) in the `dev`
  group. Each `--reinstall-package` reuses those tools from the venv
  instead of staging a fresh isolated build env — which would re-install
  ziglang (~50 MB) and friends every time, ~5s of pure overhead (~9s total
  vs ~4s).
- **no `uv cache clean`** in `reinstall`. `--reinstall-package` already
  forces a rebuild, so the clean buys nothing — and it serializes on uv's
  global cache lock, which once stalled for the full 300s lock timeout when
  another uv process held it. That 300s hang, not the build, was the only
  genuinely slow run.

`ci-test` uses `uv run --no-sync` so the already-built env (from `reinstall`
locally, or `uv sync` in CI) isn't rebuilt a second time.

An *editable* install would cut this to ~0.5s via meson-python's
rebuild-on-import, but that hook shells out to `ninja` and is fragile under
uv's build isolation (stale `ninja` path → `FileNotFoundError`); it needs
the same `no-build-isolation` setup just to work. Not worth the
rebuild-on-import machinery for a ~4s gap — see `todo.md`.

## The skar_zig dependency: local path vs. release pin

`src/zig/build.zig.zon` pins the upstream `skar` package. There are
two modes, and they behave differently for in-place builds vs.
sdist/wheel builds:

```zig
// Local development — resolve from the sibling checkout on disk.
.skar = .{ .path = "../../../skar_zig" },
```

- **In-place builds** (`just test`, `uv sync`) build from the source
  tree, so the relative `.path` resolves and everything works **with
  no network and without `skar_zig` on GitHub**. This is the loop for
  developing both sides together.
- **sdist / wheel builds** (`uv build`, `just wheel`, CI) build from
  an *isolated* sdist copy in a temp dir, where `../../../skar_zig` no
  longer exists. These require a URL+hash pin instead.

To switch to a release pin (needed before building wheels / publishing
/ running the `wheels` CI workflow):

```sh
cd src/zig && zig fetch --save=skar \
  https://github.com/ajfriend/skar_zig/archive/refs/tags/vX.Y.Z.tar.gz
```

That rewrites the `dependencies.skar` entry from `.path` to `.url` +
`.hash`. Re-run `just test` to confirm the pinned version still works.
(This is also how you bump to a newer `skar_zig` later — re-run
`zig fetch --save` against the new tag.)

## Continuous integration

- `.github/workflows/test.yml` — runs `just ci-test` across
  {Linux, macOS, Windows} × {3.11–3.14}.
- `.github/workflows/wheels.yml` — builds the sdist + the full
  cibuildwheel matrix, and on a published GitHub release publishes to
  PyPI via OIDC trusted publishing.

Both build from an isolated sdist, so **CI only works once
`build.zig.zon` uses a URL+hash pin** (see the section above) — i.e.
after `skar_zig` is pushed and tagged on GitHub. With the local
`.path` pin, the in-place `test` workflow's `uv sync` step still works,
but the `wheels` workflow's sdist build does not.

## Cutting a release

Once the `pypi` Trusted-Publisher OIDC environment is configured on the
GitHub repo (no API tokens involved):

1. **Pin a released `skar_zig`** (see above) and commit it.
2. **Bump the version** in `pyproject.toml` (`project.version`). Commit
   + push to `main`. Wait for `test` and `wheels` to go green.
3. **Create a tag + release** on GitHub: Draft a new release → choose
   tag `vX.Y.Z` ("Create new tag on publish") → target `main` →
   Generate release notes → Publish.
4. **Watch the publish.** The release event triggers the `wheels`
   workflow; its `to-pypi` job downloads every artifact and pushes to
   PyPI via OIDC. If the `pypi` environment requires reviewers, approve
   the deployment.
5. **Verify** with `pip install skar==X.Y.Z` in a fresh venv.

PyPI versions and tags are immutable — if an upload fails mid-publish,
bump to `X.Y.Z+1` rather than reusing the number.

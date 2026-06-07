# Non-editable install (the wheel is exercised as shipped). uv won't pick up a
# source change on its own, so `test` rebuilds explicitly via `reinstall`: a
# cold ReleaseFast Zig build + Cython compile. `--no-build-isolation-package` +
# the `build` group reuse the backend from the venv (~4s) instead of staging a
# fresh isolated build env each time (~9s). These are *local* flags — CI and
# `uv build` use normal isolation. No `uv cache clean` (it once stalled 300s on
# the uv lock and reuses nothing).
export UV_NO_EDITABLE := "1"
export UV_OFFLINE := "0"  # toggle on when offline to avoid failures

_:
    just --list

# Self-cleaning: drop any stale *editable* install (a lingering meson-python
# editable loader + .pth + build/ dir bakes paths into a uv build-isolation
# temp dir that gets deleted, so `import skar` later fails). Then do the cold
# non-editable reinstall. Leading `-` / `_rm` no-op when nothing's there.
reinstall:
    -uv pip uninstall skar
    just _rm '_skar_editable_loader.py'
    just _rm 'skar-editable.pth'
    just _rm build
    just _rm .zig-cache
    uv sync --reinstall-package skar --no-build-isolation-package skar --group build

# rebuild (pick up source changes) + run the suite
test: reinstall ci-test

# --no-sync: the env is already built (by `reinstall` locally, by `uv sync`
# in CI), so don't let `uv run` trigger a second redundant rebuild.
ci-test:
    uv run --no-sync pytest -q

# Build the wheel.
wheel:
    uv build

# Open JupyterLab for interactive experimentation with skar.
lab:
    uv run --group lab jupyter lab

# The DGGS scripts depend on dggal (Ecere DGGAL) for the ISEA/IVEA/rHEALPix
# systems. dggal/ecrt only publish an arch-broken macOS arm64 wheel (x86_64
# dylibs inside an arm64 wheel) and there's no native arm64 eC toolchain to
# build from source, so on Apple Silicon these run under an x86_64 (Rosetta)
# Python 3.13 where the wheels are self-consistent (Linux wheels are correct
# too). A dedicated env (.venv-dggs) keeps the native arm64 dev env untouched;
# skar is rebuilt x86_64 into it via zig.
dggs_python := "cpython-3.13-macos-x86_64"
dggs_env := ".venv-dggs"

# Build skar (x86_64) + the dggs group into the dedicated Rosetta env.
_dggs-sync:
    UV_PROJECT_ENVIRONMENT={{dggs_env}} uv sync --python {{dggs_python}} \
        --reinstall-package skar --no-build-isolation-package skar \
        --group build --group dggs

# Run the DGGS aspect-ratio survey at an H3-r9-matched resolution.
# Streams cells, solves each with skar, writes PNGs to scripts/dggs/out/.
dggs: _dggs-sync
    UV_PROJECT_ENVIRONMENT={{dggs_env}} uv run --no-sync scripts/dggs/survey.py

# Recalibrate the resolutions that match H3 r9 cell area (skar-free).
# Re-run when adding a new DGGS, then bake the result into survey.py.
calibrate: _dggs-sync
    UV_PROJECT_ENVIRONMENT={{dggs_env}} uv run --no-sync scripts/dggs/calibrate.py

# Stress test: solve millions of H3 cells across all resolutions with the
# DEFAULT solver settings and assert none return did_not_converge.
dggs-stress: _dggs-sync
    UV_PROJECT_ENVIRONMENT={{dggs_env}} uv run --no-sync scripts/dggs/dnc_stress.py

# Sweep every system across all resolutions: map the DNC boundary and flag any
# non-monotonic / unexpected did_not_converge. Writes out/dnc_sweep.png.
dnc-sweep: _dggs-sync
    UV_PROJECT_ENVIRONMENT={{dggs_env}} uv run --no-sync scripts/dggs/dnc_sweep.py

# US-state aspect ratios: geopandas -> skar.solve -> plot. Writes
# scripts/states/out/states.png.
states: reinstall
    uv run --group geo scripts/states/states.py

# Country aspect ratios: geopandas -> skar.solve -> plot. Writes
# scripts/countries/out/countries.png.
countries: reinstall
    uv run --group geo scripts/countries/countries.py

purge:
    just _rm .venv
    just _rm '*.pytest_cache'
    just _rm .DS_Store
    just _rm '*.egg-info'
    just _rm dist
    just _rm __pycache__
    just _rm uv.lock
    just _rm build
    just _rm zig-out
    just _rm .zig-cache
    just _rm .mypy_cache
    just _rm .ruff_cache
    just _rm 'libskar.*'
    uv cache clean skar

_rm pattern:
    -@find . -name "{{pattern}}" -prune -exec rm -rf {} +

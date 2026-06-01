# Non-editable install (the wheel is exercised as shipped). uv won't pick up
# a source change on its own, so `test` rebuilds explicitly via `reinstall`:
# a cold ReleaseFast Zig build + Cython compile, ~4s. The build backend is
# reused from the venv (pyproject [tool.uv] no-build-isolation) rather than
# re-staged into a fresh isolated env each time (~5s saved). No `uv cache
# clean` — it adds nothing and once stalled 300s on the uv lock.
export UV_NO_EDITABLE := "1"
export UV_OFFLINE := "0"  # toggle on when offline to avoid failures

_:
    just --list

reinstall:
    uv sync --reinstall-package skar

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

# Run the DGGS finest-resolution aspect-ratio survey (H3/S2/A5).
# Streams cells, solves each with skar, writes PNGs to scripts/dggs/out/.
dggs:
    uv run --group dggs scripts/dggs/survey.py

# US-state aspect ratios: geopandas -> skar.solve -> plot. Writes
# scripts/states/out/states.png.
states:
    uv run --group geo scripts/states/states.py

# Country aspect ratios: geopandas -> skar.solve -> plot. Writes
# scripts/countries/out/countries.png.
countries:
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

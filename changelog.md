# Changelog

Notable changes to skar. Terse by design — each entry points to the PR or
commit that carries the full detail.

## [Unreleased]

- Move the DGGS investigation (the Parquet cell-cache pipeline, survey,
  calibration, DNC checks, web viewer, and live-engine explorations) to its
  own repo, [dggs_compare](https://github.com/ajfriend/dggs_compare), which
  consumes skar as a pinned release and doubles as skar's pre-release
  regression gate (see dev.md's release runbook). skar_py keeps the small
  states/countries examples.

- Bump `skar_zig` to **v0.6.0** and expose its solver-path selection:
  `skar.solve(..., method=)` takes `'alternating'` / `'trust'` / `'auto'` and
  **defaults to `'auto'`** — upstream's alias for the recommended method,
  currently the trust-region path, which converges on the wide-angle and
  elongated inputs the alternating path structurally cannot. `Converged` /
  `DidNotConverge` gain a `.method` field recording the concrete path that
  produced the outcome. Degenerate (rank-deficient) input is now rejected as
  coplanar in preprocessing on every path, even with the near-coplanarity
  check disabled. (#16)

- Bump `skar_zig` to **v0.5.0** and adapt the C shim to its per-algorithm
  `diag` union (`outer_iters` now flows through `Diagnostics.totalIters()`).
  Solver behavior and the Python API are unchanged; upstream's experimental
  `.trust`/`.auto` solver paths are not yet exposed. (#15)
- DGGS scripts: add **IVEA7H** and refactor the DGGAL wiring to a single
  registry (`dggal_common.DGGAL_SYSTEMS` + `Adapter.iter_sample`) — a new DGGAL
  grid is now one row that calibrate/survey/dnc_sweep/validate_corners all loop
  over. IVEA7H shares ISEA7H's r10 size; slightly more circular hexagons
  (median AR ~1.15), 0 DNC across r0–r19. (#8)
- DGGS scripts: add **ISEA7H** (Ecere `dggal`/DGGAL) to the aspect-ratio survey
  + DNC sweep — the first DGGAL-backed grid, via a shared `dggal_common.py`
  adapter (median AR ~1.17 at the H3-r9-matched r10). dggal's macOS arm64 wheel
  is arch-broken (x86_64 dylibs), so the `dggs` just-targets run under an x86_64
  (Rosetta) Python 3.13 in `.venv-dggs`; native arm64 dev is untouched.
  Corners-only validated (max ΔAR ≈ 1.6e-6). (#7)
- DGGS scripts/notebook: swap pure-Python `pya5` for the Rust/PyO3 `a5_fast`
  A5 binding (~33× faster A5 cell generation; `survey.py` ~2.7s → ~0.5s).
  Verified drop-in — cell ids and boundaries match pya5 (~2e-13°). a5_fast
  has no cp314 wheel yet, so syncing the `dggs`/`lab` groups builds it from
  sdist and needs a Rust toolchain. (#5)

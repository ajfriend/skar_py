# Changelog

Notable changes to skar. Terse by design — each entry points to the PR or
commit that carries the full detail.

## [Unreleased]

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

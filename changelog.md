# Changelog

Notable changes to skar. Terse by design — each entry points to the PR or
commit that carries the full detail.

## [Unreleased]

- DGGS scripts/notebook: swap pure-Python `pya5` for the Rust/PyO3 `a5_fast`
  A5 binding (~33× faster A5 cell generation; `survey.py` ~2.7s → ~0.5s).
  Verified drop-in — cell ids and boundaries match pya5 (~2e-13°). a5_fast
  has no cp314 wheel yet, so syncing the `dggs`/`lab` groups builds it from
  sdist and needs a Rust toolchain. (#5)

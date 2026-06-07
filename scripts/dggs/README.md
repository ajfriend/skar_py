# DGGS aspect-ratio survey

For `N` random cells at a commonly-used resolution of **H3**, **S2**, and
**A5**, compute the tightest enclosing-cone aspect ratio with `skar` and
plot the results. H3 res 9 is the reference (~0.1 km², a typical working
resolution); S2 and A5 use the resolution whose cell area is closest to
an H3 r9 cell — **S2 L15** (0.76×) and **A5 r14** (1.15×).

```sh
just dggs
# or: uv run --group dggs scripts/dggs/survey.py
```

The matching resolutions are computed by `calibrate.py` (`just calibrate`)
and baked into `survey.py`'s constants; re-run it when adding a new DGGS.

For an interactive version with the same plots inline and a configurable
`N`, see [`notebooks/dggs_survey.ipynb`](../../notebooks/dggs_survey.ipynb)
(open with `just lab`).

Writes two PNGs to `scripts/dggs/out/` (gitignored) and prints a summary
stats table:

- `histograms.png` — per-system aspect-ratio distributions, shared bins.
- `extremes.png` — a 3×2 grid: each system's best (most circular) and
  worst (most elongated) cell, gnomonic-projected with its enclosing
  ellipse.

## Design

Single-pass and file-free. A generator streams one cell at a time as
`(id, (M, 3) unit-vertex array)`; each is solved immediately with
`skar.solve(...)` at the strict default `gap_tol=1e-6`, and only the
running aggregates are kept — the aspect ratios (one float per cell) plus
the two extreme cells per system. Nothing is materialized except the
final PNGs.

This is the Python port of the original Zig-repo pipeline
(`gen_cells.py → data/*.json → aspect.zig → aspect.json → plots`),
collapsed into one pass now that `skar.solve` is callable from Python.

Every cell converges at the strict `1e-6` default. (A band of H3
resolutions, r7–r10, used to stall at ~1.7e-6 and needed a relaxed
`1e-5`; skar_zig v0.2.0 fixed it — see
[`h3_gap_floor_report.md`](../../h3_gap_floor_report.md).)

Config (`N`, `SEED`, resolutions, `GAP_TOL`, …) lives in constants at the
top of `survey.py` — edit in place; no CLI args.

# DGGS aspect-ratio survey

For `N` random cells at the finest resolution of **H3**, **S2**, and
**A5**, compute the tightest enclosing-cone aspect ratio with `skar` and
plot the results.

```sh
just dggs
# or: uv run --group dggs scripts/dggs/survey.py
```

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
`skar.solve(..., gap_tol=1e-3)`, and only the running aggregates are
kept — the aspect ratios (one float per cell) plus the two extreme cells
per system. Nothing is materialized except the final PNGs.

This is the Python port of the original Zig-repo pipeline
(`gen_cells.py → data/*.json → aspect.zig → aspect.json → plots`),
collapsed into one pass now that `skar.solve` is callable from Python.

`gap_tol = 1e-3` (not skar's strict `1e-6`): at finest resolution the
sub-metre S2/A5 cells hit an f64 duality-gap floor and would otherwise
return `did_not_converge`, even though their aspect ratios are accurate.
Solving at `1e-3` keeps the distribution complete.

Config (`N`, `SEED`, resolutions, `GAP_TOL`, …) lives in constants at the
top of `survey.py` — edit in place; no CLI args.

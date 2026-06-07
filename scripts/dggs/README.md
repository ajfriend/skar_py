# DGGS aspect-ratio survey

For `N` random cells at a commonly-used resolution of **H3**, **S2**, **A5**,
and the DGGAL hex grids **ISEA7H** and **IVEA7H**, compute the tightest
enclosing-cone aspect ratio with `skar` and plot the results. H3 res 9 is the
reference (~0.1 km², a typical working resolution); the others use the
resolution whose cell area is closest to an H3 r9 cell — **S2 L15** (0.76×),
**A5 r14** (1.15×), and **ISEA7H / IVEA7H r10** (1.65× — aperture-7 steps area
by 7×, so r10 is the nearest level). ISEA7H and IVEA7H share cell sizes and
topology; they differ in projection (ISEA vs IVEA), and IVEA7H's hexagons come
out a touch more regular.

```sh
just dggs
```

The matching resolutions are computed by `calibrate.py` (`just calibrate`)
and baked into `survey.py`'s constants; re-run it when adding a new DGGS.

### Platform note (DGGAL grids)

ISEA7H and IVEA7H come from [**dggal**](https://dggal.org) (Ecere's DGGS
Abstraction Library) via `scripts/dggs/dggal_common.py`. dggal/ecrt currently
publish an **arch-broken macOS arm64 wheel** (it bundles x86_64 dylibs in an
arm64 wheel), and there's no native arm64 eC toolchain to build from source, so
the `just` DGGS targets run under an **x86_64 (Rosetta) Python 3.13** in a
separate env (`.venv-dggs`) where the wheels are self-consistent. The native
arm64 dev env is untouched. Linux wheels (x86_64 + aarch64) are correct, so
CI/Linux need no special handling. Adding a DGGAL grid is **one row** in
`dggal_common.DGGAL_SYSTEMS` (class, color, matched resolution, calibrate scan
range); the four scripts loop that registry, so the remaining grids (ISEA3H,
rHEALPix, …) are one-liners.

For an interactive version with the same plots inline and a configurable
`N`, see [`notebooks/dggs_survey.ipynb`](../../notebooks/dggs_survey.ipynb)
(open with `just lab`).

Writes two PNGs to `scripts/dggs/out/` (gitignored) and prints a summary
stats table:

- `histograms.png` — per-system aspect-ratio distributions, shared bins.
- `extremes.png` — a grid with one row per system: each system's best (most
  circular) and worst (most elongated) cell, gnomonic-projected with its
  enclosing ellipse.

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

## Convergence validation

Two scripts exhaustively check that the solver "just works" at the strict
default `gap_tol=1e-6` across every resolution of all five DGGS:

- `dnc_sweep.py` (`just dnc-sweep`) — sweeps H3 (r0–15), S2 (L0–30),
  A5 (r0–30), and ISEA7H / IVEA7H (r0–19), enumerating every cell at coarse
  resolutions and sampling heavily (S2/H3 500k, A5/ISEA7H/IVEA7H 100k per res)
  at fine ones. It locates the did-not-converge (DNC) boundary, flags any
  non-monotonic behaviour or "DNC islands", dumps offending cells to
  `out/dnc_sweep_cells.txt`, and writes `out/dnc_sweep.png` (DNC % and worst
  converged gap vs resolution).
- `dnc_stress.py` (`just dggs-stress`) — a deeper H3-only stress (every
  cell ≤ r4, all 12 pentagons per res, 500k/res above).

The survey/sweep feed skar only a cell's **corner** vertices. For the
equal-area DGGAL grids (slightly non-geodesic edges) `validate_corners.py`
checks that empirically for every grid in the registry (ISEA7H, IVEA7H):
across resolutions — including the coarsest levels and the 12 pentagons — the
aspect ratio from corners matches the ratio from edge-refined vertices to
within solver tolerance (overall max ΔAR ≈ 2e-6, at the 1e-6 gap floor). So
corners-only is exact for our purposes.

As of skar_zig **v0.4.0**, the combined coverage is **>20M cells with zero
unexpected DNCs**: H3 (r0–r15) and the DGGAL hex grids **ISEA7H / IVEA7H
(r0–r19)** stay clean throughout — their near-circular hexagons (median AR
~1.15–1.17) keep the problem well-conditioned even at the finest level — while
S2/A5 only DNC at their finest sub-metre levels (onset L28/r28), the correct,
monotonic f64 duality-gap floor, not a solver defect. (Each DGGAL grid
standalone: 1.6M cells, 0 DNC, ~2.5 min; the H3/S2/A5 portions are unchanged
and dominate the full sweep's runtime.)

# DGGS aspect-ratio explorations

Throwaway-but-kept analysis scripts that dig into the **aspect-ratio (AR)
distributions of the DGGAL grids** (ISEA7H / IVEA7H) beyond what the main
`../survey.py` reports. They grew out of one question — *why don't the DGGAL
grids have any cells near AR 1.0?* — and ended up characterizing how each
projection lays out its shape distortion (and confirming a surprising "dark
spot" feature is real, not a solver bug).

These are exploration tools, not part of the survey/sweep pipeline. They reuse
`../dggal_common.py` (the `Adapter`) and `skar`, and like the rest of `dggs/`
they run under the **x86_64 (Rosetta) Python 3.13** env — see the "Platform
note" in [`../README.md`](../README.md). Each writes PNGs to
`explorations/out/` (gitignored). No CLI args — edit the in-file constants.

Run any of them with:

```sh
UV_PROJECT_ENVIRONMENT=.venv-dggs uv run --no-sync \
    scripts/dggs/explorations/<script>.py
```

(Or `just _dggs-sync` once first to (re)build the env.)

## The scripts

| script | what it shows |
|--------|---------------|
| `ar_histograms.py` | Large-sample (N points) per-grid AR histograms at the matched resolution. ISEA7H is **right-skewed** (long tail to ~1.36); IVEA7H is a **flat-topped, sharply-bounded band** (~1.06–1.22). Bump `N` for a crisper picture. |
| `ar_heatmap.py` | The AR field mapped spatially. **Global** equirectangular (the 20-fold icosahedral pattern: low-distortion face centroids + bright high-distortion seams for ISEA; a near-uniform field for IVEA) and a **single-face** gnomonic zoom. |
| `dark_spots_locate.py` | Finds the rare **sharp low-AR "dark spots"** uniform sampling misses, confirms every solve converges, and zooms one to show it's razor-sharp (AR ~1.0 at a point, ~1.34 a hair away), sitting *on* the seams. |
| `dark_spots_verify.py` | Cross-checks skar's AR at a dark spot against an **independent minimum-enclosing-ellipse** (Khachiyan MVEE). They agree to ~5 decimals → **not a skar bug**: the cell is an irregular hexagon whose bounding ellipse just happens to be circular. |
| `dark_spots_cells.py` | Plots a dark-spot cell + neighbors (orthographic) with each cell's enclosing ellipse. The grid is a valid continuous tiling of irregular hexagons; the spike cell's ellipse is circular while its neighbors' are all elongated along the seam axis. |

## What we learned

- For a small cell, the enclosing-cone AR ≈ the projection's **local Tissot
  anisotropy ratio** (the min-enclosing ellipse is affine-equivariant), so it's
  resolution-independent and the histograms/heatmaps are really pictures of each
  projection's shape-distortion field.
- **ISEA7H** (Snyder equal-area) concentrates distortion into sharp piecewise
  **seams** (sub-triangle boundaries) → a long high tail and the dark-spot
  artifacts. **IVEA7H** spreads it smoothly and **caps it** (max ~1.22) → the
  flat band. By this metric IVEA7H is the better-conditioned grid.
- Both grids *do* reach near-AR-1.0 at their ~20 face centroids (≈1.0012), but
  those are measure-≈-zero by area, so uniform sampling shows them only as a
  thin tail (or misses them).
- The "dark spots" embedded in the ISEA seams are **real geometry, not a skar
  bug** — confirmed against an independent MVEE.

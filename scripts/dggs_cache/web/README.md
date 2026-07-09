# DGGS aspect-ratio web explorer

An interactive view of the aspect-ratio survey: dynamic, overlaid histograms for
any system/resolution, and two synced orthographic globes with cells colored by
aspect ratio. A static page — no build step, no server framework. The globes are
rendered by [ajglobe](https://github.com/ajfriend/ajglobe) (vendored under
`vendor/`; refresh with `just web-vendor`).

```sh
just web        # build the data (if needed) + serve at http://localhost:8000
just web-data   # just (re)build the data
```

`just web-data` runs `gen-cells` first if you haven't (the Parquet cache).

## How it works

`build_data.py` is a third reader of the Parquet cell cache (alongside
`../survey.py` / `../dnc_check.py`): it solves every cached cell with `skar`
(native, no DGGS library) and emits browser-friendly data into `out/`
(gitignored via `scripts/**/out/`):

- `histograms.json` — for every `(system, resolution)`, a **fixed** fine-bin
  histogram (`NBINS` bins over `[1, amax]` + an overflow bin) plus summary
  stats. The fixed grid lets the page re-aggregate to any coarser bin width in
  the browser; `amax` is a high global percentile so the working resolutions
  never reach overflow.
- `globe/{sys}_r{res}_{pos.f32,idx.u32,ar.f32,ids.json}` — the coarse
  resolutions only (largest res per system with ≤ `GLOBE_MAX_CELLS` cells, and
  everything below). ajglobe's native flat-binary polygon format: Float32
  `[lng, lat]` vertex pairs, Uint32 ring starts, Float32 AR per cell
  (NaN = did-not-converge), and the cell-id strings for hover. Open rings, any
  winding — ajglobe triangulates by ring topology, so there's no orientation or
  antimeridian preprocessing.
- `manifest.json` — what exists, per-system web colors (matched to the survey
  PNGs), the bin grid, and the shared globe AR max.

The page (`index.html` + `app.js`; d3 + Observable Plot from a CDN for the
histograms, vendored ajglobe for the globes):

- **Histograms** (Observable Plot) — pick any number of `system · resolution`
  series; one shared aspect-ratio axis, a bins slider, density/count, and
  linear/log. The displayed x-domain tracks the current selection.
- **Globe** (ajglobe `Orb`, WebGL2) — two GPU-rendered globes that rotate and
  zoom together (drag either; arrows/WASD with a globe focused), cells colored
  by AR, GPU hover picking (cell id + AR), Natural Earth country outlines for
  context. Color domain is **per-globe** by default (each globe's own AR
  quantiles, to show within-system structure) or **shared** (one absolute
  scale, to compare systems); toggling restyles in place — no re-tessellation.

Edit the constants at the top of `build_data.py` (no CLI args, project
convention) to change `NBINS`, the `amax` percentile, or the globe cell cap.

## Full-globe experiment

`globe_full.html` renders **every** ivea7h cell (up to 1.18M at r6) rather than
the sampled cache — same binary format, produced by the two-pass pipeline
(`just web-full-geom` for DGGAL geometry under Rosetta, then `just web-full`
for native skar ARs) into `out/full/`. Both passes skip levels already on disk.

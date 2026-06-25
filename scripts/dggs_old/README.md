# DGGS live-engine analyses

The DGGS analyses that **can't** run off the Parquet cell cache — they query the
live DGGS engine for things a bag-of-cells snapshot doesn't capture (arbitrary
point→cell over a grid, a cell's neighbors, edge-refined boundaries). Everything
that *can* read the cache — the survey, calibrate, the DNC tests, and the AR
distribution histograms — moved to [`../dggs_cache/`](../dggs_cache/).

These run under the x86_64 (Rosetta) `.venv-dggs` env, because dggal ships an
arch-broken macOS arm64 wheel (see `dggal_common.py` and the platform note in
`expansion_plan.md`):

```sh
UV_PROJECT_ENVIRONMENT=.venv-dggs uv run --no-sync \
    scripts/dggs_old/explorations/ar_heatmap.py
```

## Contents

- **`dggal_common.py`** — the live DGGAL engine wrapper: `Adapter` (point→cell,
  vertices, centroids, neighbors, edge refinement), the `DGGAL_SYSTEMS` registry,
  and `latlng_ring`. (`dggs_cache` carries its own copy, used only by `gen_dggal`
  to *produce* the Parquet.)
- **`validate_corners.py`** — confirms the corners-only AR metric is exact:
  compares each cell's corner-derived ratio against an edge-refined reference
  (`getZoneRefinedWGS84Vertices`), per resolution, for ISEA7H/IVEA7H.
- **`explorations/`** — `ar_heatmap.py` (AR over a lon/lat grid + per-face
  gnomonic zoom), `dark_spots_locate.py` (grid scan for rare low-AR seam cells),
  `dark_spots_verify.py` (perturb a point and re-query), `dark_spots_cells.py`
  (a spike cell + its neighbors).
- **`expansion_plan.md`** — adding more DGGAL grids (one `DGGAL_SYSTEMS` row).

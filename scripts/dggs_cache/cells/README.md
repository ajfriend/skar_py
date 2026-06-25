# Random cell generators

One standalone script per DGGS that generates random cells and writes them to
Parquet files in a common schema. Generation (which needs the per-DGGS library)
is decoupled from analysis (which needs `skar`): an analysis reads the Parquet
back with numpy + pyarrow only, so it never has to import — or satisfy the
version/arch constraints of — h3, s2sphere, a5_fast, or dggal.

Each generator writes **one file per resolution**, from `0` up to the system's
finest. The per-resolution budget is dense up to the system's working (target)
resolution — `N_BIG` (~100k) cells, where the survey and AR explorations want
lots of cells — and thin beyond it — `N_SMALL` (~25k), where only the
`dnc_sweep`/`calibrate` tests read and coverage matters more than count.
(Keeping the deep tail thin is what keeps generation fast; dggal samples every
cell under Rosetta.) At each resolution, if the whole resolution has `<= n`
cells it enumerates them all (exact, complete — coarse resolutions saturate long
before `n` random samples would); otherwise it draws `n` uniform-on-sphere
points and dedups. So every consumer just reads whatever cells exist at the
resolution it wants.

Each generator is a PEP 723 / `uv run` script carrying its own dependency and
Python version, so the libraries never have to coexist in one environment. The
files are cacheable: each `(dggs, res)` maps to one file
`out/{dggs}_r{res}.parquet` (gitignored), e.g. `h3_r9.parquet`. The budgets
(`N_BIG`/`N_SMALL`) and `SEED` are pipeline config in `_common.py`, not encoded
in the filename. A file holds `min(n, distinct samples)` cells, or all `count`
cells when the resolution was enumerated.

## Schema

One row per distinct cell:

| column  | type                              | notes                          |
| ------- | --------------------------------- | ------------------------------ |
| `dggs`  | `string`                          | constant per file, e.g. `h3`   |
| `res`   | `int32`                           | constant per file (resolution) |
| `cid`   | `string`                          | cell id text                   |
| `verts` | `list<fixed_size_list<double,2>>` | ring of `[lat, lng]` degrees   |

`dggs`/`res` are constant per file (a few bytes after RLE) but kept as columns
so files concatenate into one self-describing table you can `groupby`. `verts`
is an open ring (no repeated closing vertex); length varies by cell (6 for H3
hexagons, 5 for A5 pentagons, 4 for S2 quads, plus a few extra at icosahedron
edges).

Rows are sorted by `cid` for a canonical, deterministic order (independent of
sampling order) — which also enables Parquet `cid` page-stats / range pushdown.
Files are written with zstd, `BYTE_STREAM_SPLIT` on the vertex floats, and
`DELTA_BYTE_ARRAY` on the sorted `cid`s — a lossless ~30% shrink over the
snappy/plain default (the lat/lng values share sign/exponent bytes that
compress; coordinates stay exact `float64`). Any modern Parquet reader
(pyarrow, DuckDB, polars) decodes this transparently.

## Run

One `uv run` per system writes all its resolutions:

```sh
uv run scripts/dggs_cache/cells/gen_h3.py     # -> h3_r0.parquet .. h3_r15.parquet
uv run scripts/dggs_cache/cells/gen_s2.py
uv run scripts/dggs_cache/cells/gen_a5.py
uv run scripts/dggs_cache/cells/gen_dggal.py  # isea7h + ivea7h
```

DGGAL ships an arch-broken macOS arm64 wheel, so on Apple Silicon `gen_dggal`
re-execs itself under an x86_64 (Rosetta) Python 3.13 (where the wheel is
self-consistent); the command above is the same on every platform (Linux wheels
are correct, so the re-exec is a no-op there).

Each generator's `MAX_RES`/`MAX_LEVEL` is at the top of its script; `N_BIG`,
`N_SMALL`, `SEED`, and the per-system `TARGET_RES` are pipeline config in
`_common.py` — no CLI args, per project convention.

## Read (analysis side)

```python
import _common  # or replicate cells_path()/the read

for cid, verts in _common.load_cells('h3', 9):
    # verts: (M, 2) array of [lat, lng] degrees
    r = skar.solve(skar.to_vec3(verts, geo='latlng_deg'), geo='vec3')
```

`_common.py` holds the shared sampler, schema, path convention, and read/write
helpers; the generators only map `(lng, lat)` to a cell id and its vertex ring.

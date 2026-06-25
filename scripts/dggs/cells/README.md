# Random cell generators

One standalone script per DGGS that generates random cells and writes them to
Parquet files in a common schema. Generation (which needs the per-DGGS library)
is decoupled from analysis (which needs `skar`): an analysis reads the Parquet
back with numpy + pyarrow only, so it never has to import — or satisfy the
version/arch constraints of — h3, s2sphere, a5_fast, or dggal.

Each generator writes **one file per resolution**, `0` up to its target (the
resolution near H3 r9's cell size). At each resolution it picks the right
strategy automatically: if the whole resolution has `<= N` cells it enumerates
them all (exact, complete — coarse resolutions saturate long before `N` random
samples would, and sampling would miss the tail); otherwise it draws `N`
uniform-on-sphere points and dedups to distinct cells.

Each generator is a PEP 723 / `uv run` script carrying its own dependency and
Python version, so the libraries never have to coexist in one environment. The
files are cacheable: a given `(dggs, res, n, seed)` maps to one file under
`out/` (gitignored); reuse it across analyses, or delete it to regenerate. `n`
is the sample budget — a file holds `min(count, distinct samples)` cells, or all
`count` cells when the resolution was enumerated.

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

```sh
uv run scripts/dggs/cells/gen_h3.py     # -> out/h3_r{0..9}_n100000_s00c0ffee.parquet
uv run scripts/dggs/cells/gen_s2.py
uv run scripts/dggs/cells/gen_a5.py
```

DGGAL ships an arch-broken macOS arm64 wheel, so on Apple Silicon run its
generator under an x86_64 (Rosetta) Python 3.13 (Linux wheels are fine with a
plain `uv run`):

```sh
uv run --python cpython-3.13-macos-x86_64 scripts/dggs/cells/gen_dggal.py
```

Knobs (`TARGET_RES`/`TARGET_LEVEL`, `N`, `SEED`) are edited in place at the top
of each script — no CLI args, per project convention.

## Read (analysis side)

```python
import _common  # or replicate cells_path()/the read

for cid, verts in _common.load_cells('h3', 9, 100_000, 0xC0FFEE):
    # verts: (M, 2) array of [lat, lng] degrees
    r = skar.solve(skar.to_vec3(verts, geo='latlng_deg'), geo='vec3')
```

`_common.py` holds the shared sampler, schema, path convention, and read/write
helpers; the generators only map `(lng, lat)` to a cell id and its vertex ring.

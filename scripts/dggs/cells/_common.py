"""Shared glue for the per-DGGS cell generators in this folder.

Each generator (``gen_<dggs>.py``) is a standalone PEP 723 / ``uv run`` script
that carries its *own* DGGS-library dependency (h3, s2sphere, a5_fast, dggal)
and Python version. This module holds the parts they share — uniform-on-sphere
sampling, the Parquet schema + writer, and the on-disk path convention — so a
generator only has to map ``(lng, lat)`` points to a cell id and that cell's
vertex ring.

The generators are deliberately ``skar``-free: they emit raw geometry, not
solved aspect ratios. Analyses read the Parquet back (numpy + pyarrow only,
native arch) and do the solving, fully decoupled from the DGGS libraries.

Schema (one row per distinct cell):

    dggs   string                            constant per file, e.g. 'isea7h'
    res    int32                             constant per file, cell resolution/level
    cid    string                            cell id text
    verts  list<fixed_size_list<double, 2>>  ring of [lat, lng] degrees (open)

Files land in ``out/`` next to this module (gitignored via ``scripts/**/out/``),
named ``{dggs}_r{res}_{kind}.parquet`` where ``kind`` is ``big`` or ``small``
(the two density tiers — see ``generate_big_small``). Read back with
``load_cells(dggs, res, kind)`` or any Parquet reader (pandas, DuckDB). SEED and
the per-kind N are pipeline config below, not encoded in the filename.
"""

import os
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

OUT_DIR = Path(__file__).resolve().parent / 'out'

# ----- pipeline config (single source of truth) --------------------------
# These pin the cache keys, so every generator and analysis must agree on them
# — hence they live here, not restated per script. SEED + N must match between
# a generator and the analysis that reads its files.
SEED = 0xC0FFEE
N_BIG, N_SMALL = 100_000, 25_000     # cells/resolution in the big / small sets
N_BY_KIND = {'big': N_BIG, 'small': N_SMALL}   # filename tag -> sample budget
# Working ("target") resolution per system: the finest in actual use, matched
# to an H3 r9 cell by calibrate.py. The big set spans 0..target; survey and
# dnc_stress read it. Update here when calibrate picks a new value.
TARGET_RES = {'h3': 9, 's2': 15, 'a5': 14, 'isea7h': 10, 'ivea7h': 10}
# -------------------------------------------------------------------------

# fixed_size_list(2): each vertex is exactly [lat, lng] deg; the outer list is
# the variable-length ring (6 for hexagons, 5 for the pentagons, 4 for quads).
VERTS_TYPE = pa.list_(pa.list_(pa.float64(), 2))
# Leaf column path of the nested verts doubles (for per-column encoding).
VERTS_LEAF = 'verts.list.element.list.element'
SCHEMA = pa.schema([
    ('dggs', pa.string()),
    ('res', pa.int32()),
    ('cid', pa.string()),
    ('verts', VERTS_TYPE),
])


def sample_uniform_lnglat(n, rng):
    """Uniform-on-sphere samples as (lng_deg, lat_deg), shape (n, 2)."""
    lng = 360.0 * rng.random(n) - 180.0
    lat = np.degrees(np.arcsin(2.0 * rng.random(n) - 1.0))  # equal-area in lat
    return np.column_stack([lng, lat])


def open_ring(ring):
    """Drop a closing vertex if the ring repeats its first point (corners only)."""
    if len(ring) >= 2 and tuple(ring[0]) == tuple(ring[-1]):
        return ring[:-1]
    return ring


def cells_path(dggs, res, kind):
    """Canonical Parquet path for a generated cell set (kind: 'big' | 'small')."""
    return OUT_DIR / f'{dggs}_r{res}_{kind}.parquet'


def generate(dggs, res, kind, *, latlng_to_cell, cid_str, cell_boundary,
             count_at=None, enumerate_at=None):
    """Build one `(dggs, res)` cell set of the given `kind` and write it.

    Draws `N_BY_KIND[kind]` cells (with the module SEED). If the resolution is
    small enough — `enumerate_at` given and `count_at(res) <= N` — every cell is
    enumerated (exact, complete; the coarse resolutions saturate well before `N`
    random samples would, and pure sampling would miss the tail). Otherwise `N`
    uniform-on-sphere points are drawn and deduped to distinct cells.

    Callbacks (one DGGS library each):
        latlng_to_cell(res, lat, lng) -> z   native, hashable cell id for the point
        cid_str(z)             -> str text id stored in the file
        cell_boundary(z)             -> sequence of (lat, lng) degree pairs (open ring)
        count_at(res)          -> int total cells at the resolution (optional)
        enumerate_at(res)      -> iterable of native ids (optional)

    Returns the written path.
    """
    n = N_BY_KIND[kind]
    if enumerate_at is not None and count_at is not None and count_at(res) <= n:
        zones = list(enumerate_at(res))
        mode = 'all'
    else:
        rng = np.random.default_rng(SEED)
        seen, zones = set(), []
        for lng, lat in sample_uniform_lnglat(n, rng):
            z = latlng_to_cell(res, float(lat), float(lng))
            if z not in seen:
                seen.add(z)
                zones.append(z)
        mode = 'sample'

    # Sort by cid for a canonical, deterministic row order (independent of
    # sampling order): enables Parquet cid page-stats / range pushdown, and lets
    # DELTA_BYTE_ARRAY prefix-compress the now-sorted ids. cell_boundary is called
    # once per cell here.
    rows = sorted(((cid_str(z), open_ring(cell_boundary(z))) for z in zones),
                  key=lambda cr: cr[0])
    cids = [c for c, _ in rows]
    verts = [[[float(la), float(ln)] for la, ln in ring] for _, ring in rows]
    table = pa.table({
        'dggs': pa.array([dggs] * len(cids), pa.string()),
        'res': pa.array([res] * len(cids), pa.int32()),
        'cid': pa.array(cids, pa.string()),
        'verts': pa.array(verts, VERTS_TYPE),
    }, schema=SCHEMA)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = cells_path(dggs, res, kind)
    # BYTE_STREAM_SPLIT packs the vertex float64s ~28% smaller, losslessly: the
    # lat/lng values share sign/exponent bytes that zstd then compresses, while
    # the random mantissa bytes are kept out of the way. DELTA_BYTE_ARRAY
    # prefix-compresses the sorted cids; the constant columns keep dictionary.
    pq.write_table(
        table, path,
        compression='zstd',
        use_dictionary=['dggs', 'res'],
        column_encoding={'cid': 'DELTA_BYTE_ARRAY', VERTS_LEAF: 'BYTE_STREAM_SPLIT'},
    )
    kb = os.path.getsize(path) / 1024
    print(f'[{dggs} r{res:<2}] {mode:>6} {len(cids):>7} cells '
          f'-> {path.name} ({kb:.0f} KiB)')
    return path


def generate_levels(dggs, max_res, kind, **callbacks):
    """Write one Parquet cell set per resolution 0..max_res (inclusive)."""
    for res in range(max_res + 1):
        generate(dggs, res, kind, **callbacks)


def generate_big_small(dggs, target_res, max_res, **callbacks):
    """Write a system's two cell sets:

    - "big":   resolutions 0..target_res at N_BIG — the dense working-resolution
               set the survey and AR explorations read.
    - "small": resolutions 0..max_res at N_SMALL — a thin all-resolution set for
               calibrate (area scan) and the DNC sweep/stress tests, which need
               every resolution (incl. the finest) but tolerate small N.
    """
    generate_levels(dggs, target_res, 'big', **callbacks)
    generate_levels(dggs, max_res, 'small', **callbacks)


def available_resolutions(dggs, kind):
    """Sorted resolutions that have a cached `(dggs, kind)` cell set."""
    import re
    suffix = f'_{kind}.parquet'
    pat = re.compile(rf'^{re.escape(dggs)}_r(\d+){re.escape(suffix)}$')
    res = [int(m.group(1)) for p in OUT_DIR.glob(f'{dggs}_r*{suffix}')
           if (m := pat.match(p.name))]
    return sorted(res)


def load_cells(dggs, res, kind):
    """Yield (cid, (M, 2) lat/lng array) for a generated cell set."""
    path = cells_path(dggs, res, kind)
    if not path.exists():
        raise FileNotFoundError(
            f'{path} not found — generate the cell sets first with '
            f'`just gen-cells` (or the matching scripts/dggs/cells/gen_*.py).')
    # Stream in batches (only the columns we need) so memory stays flat — the
    # analyses solve one cell at a time and never need the whole table at once.
    for batch in pq.ParquetFile(path).iter_batches(columns=['cid', 'verts']):
        for row in batch.to_pylist():
            yield row['cid'], np.asarray(row['verts'], dtype=float)

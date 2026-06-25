# /// script
# requires-python = ">=3.11"
# dependencies = ["a5_fast", "numpy>=1.24", "pyarrow>=15"]
# ///
"""Generate random A5 cell sets, one per resolution 0..TARGET_RES -> Parquet.

Coarse resolutions (<= N cells) are enumerated in full; finer ones are sampled.
Schema + write logic in _common.py.

Run:  uv run scripts/dggs/cells/gen_a5.py
No CLI args (project convention) — edit the constants below in place.
"""

import a5_fast as a5  # Rust/PyO3 A5 binding (~30x faster than pure-Python pya5)

import _common

# ----- knobs -------------------------------------------------------------
TARGET_RES = 14         # a5 supports 0..30 (a5.MAX_RESOLUTION); r14 ~0.13 km^2
N = 100_000
SEED = 0xC0FFEE
# -------------------------------------------------------------------------


def zone_at(res, lon, lat):
    return a5.lonlat_to_cell(lon, lat, res)  # int, hashable


def count_at(res):
    return a5.get_num_cells(res)


def enumerate_at(res):
    res0 = a5.get_res0_cells()
    if res == 0:
        yield from res0
    else:
        for c0 in res0:
            yield from a5.cell_to_children(c0, res)


def ring_of(cid):
    ring = a5.cell_to_boundary(cid)  # closed ring of (lon, lat)
    return [(lat, lon) for lon, lat in ring]  # -> (lat, lon); closing repeat dropped upstream


if __name__ == '__main__':
    _common.generate_levels(
        'a5', TARGET_RES, N, SEED,
        zone_at=zone_at, cid_str=a5.u64_to_hex, ring_of=ring_of,
        count_at=count_at, enumerate_at=enumerate_at)

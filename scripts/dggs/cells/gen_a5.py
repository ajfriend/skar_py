# /// script
# requires-python = ">=3.11"
# dependencies = ["a5_fast", "numpy>=1.24", "pyarrow>=15"]
# ///
"""Generate random A5 cell sets -> Parquet, one file per resolution 0..MAX_RES.

Coarse resolutions (<= N cells) are enumerated in full; finer ones sample N.
Schema + write logic (and N/SEED config) in _common.py.

Run:  uv run scripts/dggs/cells/gen_a5.py
No CLI args (project convention) — edit the constants below in place.
"""

import a5_fast as a5  # Rust/PyO3 A5 binding (~30x faster than pure-Python pya5)

import _common

# ----- knobs -------------------------------------------------------------
MAX_RES = 30            # finest a5 resolution to generate (a5.MAX_RESOLUTION)
# -------------------------------------------------------------------------


def latlng_to_cell(res, lat, lng):
    return a5.lonlat_to_cell(lng, lat, res)  # int, hashable


def count_at(res):
    return a5.get_num_cells(res)


def enumerate_at(res):
    res0 = a5.get_res0_cells()
    if res == 0:
        yield from res0
    else:
        for c0 in res0:
            yield from a5.cell_to_children(c0, res)


def cell_boundary(cid):
    ring = a5.cell_to_boundary(cid)  # closed ring of (lng, lat)
    return [(lat, lng) for lng, lat in ring]  # -> (lat, lng); closing repeat dropped upstream


if __name__ == '__main__':
    _common.generate_levels(
        'a5', MAX_RES,
        latlng_to_cell=latlng_to_cell, cid_str=a5.u64_to_hex, cell_boundary=cell_boundary,
        count_at=count_at, enumerate_at=enumerate_at)

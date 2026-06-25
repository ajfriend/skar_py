# /// script
# requires-python = ">=3.11"
# dependencies = ["h3>=4", "numpy>=1.24", "pyarrow>=15"]
# ///
"""Generate random H3 cell sets -> Parquet, one file per resolution 0..MAX_RES.

Coarse resolutions (<= N cells) are enumerated in full; finer ones sample N.
Schema + write logic (and N/SEED config) in _common.py.

Run:  uv run scripts/dggs/cells/gen_h3.py
No CLI args (project convention) — edit the constants below in place.
"""

import h3

import _common

# ----- knobs -------------------------------------------------------------
MAX_RES = 15            # finest h3 resolution to generate (h3 supports 0..15)
# -------------------------------------------------------------------------


def latlng_to_cell(res, lat, lng):
    return h3.latlng_to_cell(lat, lng, res)


def count_at(res):
    return h3.get_num_cells(res)


def enumerate_at(res):
    res0 = h3.get_res0_cells()
    if res == 0:
        yield from res0
    else:
        for c0 in res0:
            yield from h3.cell_to_children(c0, res)


def cell_boundary(cid):
    return h3.cell_to_boundary(cid)  # [(lat, lng), ...] degrees, corners only


if __name__ == '__main__':
    _common.generate_levels(
        'h3', MAX_RES,
        latlng_to_cell=latlng_to_cell, cid_str=str, cell_boundary=cell_boundary,
        count_at=count_at, enumerate_at=enumerate_at)

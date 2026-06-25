# /// script
# requires-python = ">=3.11"
# dependencies = ["h3>=4", "numpy>=1.24", "pyarrow>=15"]
# ///
"""Generate random H3 cell sets -> Parquet: a dense "big" set at the working
resolutions and a thin "small" set at every resolution.

Coarse resolutions (<= N cells) are enumerated in full; finer ones are sampled.
Schema + write logic in _common.py.

Run:  uv run scripts/dggs/cells/gen_h3.py
No CLI args (project convention) — edit the constants below in place.
"""

import h3

import _common

# ----- knobs -------------------------------------------------------------
TARGET_RES = 9          # h3 supports 0..15; r9 ~0.1 km^2 (survey reference)
MAX_RES = 15            # finest h3 resolution (for the all-res small set)
N_BIG = 100_000         # dense set, 0..TARGET_RES (survey, AR explorations)
N_SMALL = 25_000        # thin all-res set, 0..MAX_RES (calibrate, DNC tests)
SEED = 0xC0FFEE
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
    _common.generate_big_small(
        'h3', TARGET_RES, MAX_RES, N_BIG, N_SMALL, SEED,
        latlng_to_cell=latlng_to_cell, cid_str=str, cell_boundary=cell_boundary,
        count_at=count_at, enumerate_at=enumerate_at)

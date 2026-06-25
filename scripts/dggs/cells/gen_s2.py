# /// script
# requires-python = ">=3.11"
# dependencies = ["s2sphere", "numpy>=1.24", "pyarrow>=15"]
# ///
"""Generate random S2 cell sets, one per level 0..TARGET_LEVEL -> Parquet.

Coarse levels (<= N cells) are enumerated in full; finer ones are sampled.
Schema + write logic in _common.py.

Run:  uv run scripts/dggs/cells/gen_s2.py
No CLI args (project convention) — edit the constants below in place.
"""

import s2sphere

import _common

# ----- knobs -------------------------------------------------------------
TARGET_LEVEL = 15       # s2sphere supports 0..30; L15 ~0.083 km^2 (survey match)
N = 100_000
SEED = 0xC0FFEE
# -------------------------------------------------------------------------


def latlng_to_cell(res, lat, lng):
    cid = s2sphere.CellId.from_lat_lng(s2sphere.LatLng.from_degrees(lat, lng))
    if res != 30:
        cid = cid.parent(res)
    return cid.id()  # int, hashable


def count_at(res):
    return 6 * 4 ** res


def enumerate_at(res):
    for cid in s2sphere.CellId.walk(res):
        yield cid.id()


def cid_str(zid):
    return format(zid, '016x')


def cell_boundary(zid):
    cell = s2sphere.Cell(s2sphere.CellId(zid))
    ring = []
    for i in range(4):
        ll = s2sphere.LatLng.from_point(cell.get_vertex(i))
        ring.append((ll.lat().degrees, ll.lng().degrees))
    return ring


if __name__ == '__main__':
    _common.generate_levels(
        's2', TARGET_LEVEL, N, SEED,
        latlng_to_cell=latlng_to_cell, cid_str=cid_str, cell_boundary=cell_boundary,
        count_at=count_at, enumerate_at=enumerate_at)

# /// script
# requires-python = ">=3.11"
# dependencies = ["s2sphere", "numpy>=1.24", "pyarrow>=15"]
# ///
"""Generate random S2 cell sets -> Parquet, one file per level 0..MAX_LEVEL.

Coarse levels (<= N cells) are enumerated in full; finer ones sample N.
Schema + write logic (and N/SEED config) in _common.py.

Run:  uv run scripts/dggs_cache/cells/gen_s2.py
No CLI args (project convention) — edit the constants below in place.
"""

import s2sphere

import _common

# ----- knobs -------------------------------------------------------------
MAX_LEVEL = 30          # finest s2 level to generate (s2sphere supports 0..30)
# -------------------------------------------------------------------------


def latlng_to_cell(res, lat, lng):
    # from_lat_lng yields a level-30 leaf; walk up to the cell at `res`.
    leaf = s2sphere.CellId.from_lat_lng(s2sphere.LatLng.from_degrees(lat, lng))
    return leaf.parent(res).id()  # int, hashable


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
        's2', MAX_LEVEL,
        latlng_to_cell=latlng_to_cell, cid_str=cid_str, cell_boundary=cell_boundary,
        count_at=count_at, enumerate_at=enumerate_at)

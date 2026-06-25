"""Match S2/A5 resolutions to an H3-res-9 cell by area.

Picks, for each DGGS, the resolution whose median cell area is closest
(in log-ratio) to the reference H3 res-9 cell. Used to choose the
resolutions baked into survey.py so all three systems compare cells of
roughly the same size.

Skar-free by design: areas come from `sparea` (spherical-polygon area),
not from the solver, so this runs even if the skar build is broken. The
random sampling mirrors survey.py (same uniform-on-sphere sampler and
SEED) so it measures the same cells the survey solves.

Adding a new DGGS = add one `<sys>_area(res, n)` function returning a
median cell area in km^2, register it in AREA_FN, and give it a SCAN
range. Then re-run and bake the printed pick into survey.py.

Run with:  just calibrate   (or: uv run --group dggs scripts/dggs/calibrate.py)
No CLI args (project convention) — edit the constants below in place.
"""

from functools import partial

import numpy as np

import a5_fast as a5  # Rust/PyO3 A5 binding (~30x faster than pure-Python pya5)
import h3
import s2sphere

import sparea

import dggal_common  # Ecere DGGAL binding glue (ISEA/IVEA hex, rHEALPix, ...)

# ----- knobs -------------------------------------------------------------
N = 2_000                   # cells sampled per resolution (median is robust)
SEED = 0xC0FFEE             # match survey.py for reproducibility
R_KM = 6371.0088            # mean Earth radius; steradian -> km^2 is R^2
SR2KM2 = R_KM * R_KM

TARGET = ('h3', 9)          # reference system + resolution
SCAN = {'s2': range(10, 20), 'a5': range(8, 20)}  # candidate resolutions
# (DGGAL systems add their own SCAN ranges from the registry below.)
# -------------------------------------------------------------------------


# ----- per-system median cell area (km^2) over N random cells ------------
# Each samples from the same SEED so candidate resolutions are comparable.
def h3_area(res, n):
    rng = np.random.default_rng(SEED)
    a = []
    for lon, lat in dggal_common.sample_uniform_lnglat(n, rng):
        cid = h3.latlng_to_cell(float(lat), float(lon), res)
        b = h3.cell_to_boundary(cid)  # [(lat, lng), ...] deg
        a.append(sparea.area(b, geo='latlng'))
    return float(np.median(a)) * SR2KM2


def s2_area(res, n):
    rng = np.random.default_rng(SEED)
    a = []
    for lon, lat in dggal_common.sample_uniform_lnglat(n, rng):
        cid = s2sphere.CellId.from_lat_lng(
            s2sphere.LatLng.from_degrees(float(lat), float(lon))).parent(res)
        cell = s2sphere.Cell(cid)
        v = np.array([tuple(cell.get_vertex(i))[:3] for i in range(4)], dtype=float)
        v /= np.linalg.norm(v, axis=1, keepdims=True)
        a.append(sparea.area(v, geo='vec3'))
    return float(np.median(a)) * SR2KM2


def a5_area(res, n):
    rng = np.random.default_rng(SEED)
    a = []
    for lon, lat in dggal_common.sample_uniform_lnglat(n, rng):
        cid = a5.lonlat_to_cell(float(lon), float(lat), res)
        ring = a5.cell_to_boundary(cid)  # closed ring of (lon, lat)
        if len(ring) >= 2 and tuple(ring[0]) == tuple(ring[-1]):
            ring = ring[:-1]
        latlng = [(lat_, lon_) for lon_, lat_ in ring]
        a.append(sparea.area(latlng, geo='latlng'))
    return float(np.median(a)) * SR2KM2


AREA_FN = {'h3': h3_area, 's2': s2_area, 'a5': a5_area}


# DGGAL systems: register an area fn + scan range from each registry row, so
# adding a grid is one line in dggal_common.DGGAL_SYSTEMS (no per-system code).
for _k, _s in dggal_common.DGGAL_SYSTEMS.items():
    AREA_FN[_k] = partial(dggal_common.Adapter(_s['cls']).area_km2, seed=SEED)
    SCAN[_k] = _s['scan']


def main():
    tsys, tres = TARGET
    target = AREA_FN[tsys](tres, N)
    print(f'target: {tsys} r{tres} median area = {target:.6f} km^2  '
          f'(N={N}, seed={SEED:#x})\n')

    for sys, scan in SCAN.items():
        rows = [(res, AREA_FN[sys](res, N)) for res in scan]
        best = min(rows, key=lambda r: abs(np.log(r[1] / target)))
        print(f'--- {sys} (target {tsys} r{tres}) ---')
        print(f'{"res":>4} {"area_km2":>12} {"ratio":>8}')
        for res, area in rows:
            mark = '  <== pick' if res == best[0] else ''
            print(f'{res:>4} {area:>12.6f} {area / target:>8.3f}{mark}')
        print(f'-> {sys} r{best[0]}  ({best[1] / target:.3f}x {tsys} r{tres})\n')


if __name__ == '__main__':
    main()

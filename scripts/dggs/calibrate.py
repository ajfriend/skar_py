# /// script
# requires-python = ">=3.11"
# dependencies = ["sparea>=0.4", "numpy>=1.24", "pyarrow>=15"]
# ///
"""Match S2 / A5 / DGGAL resolutions to an H3-res-9 cell by area.

Picks, for each DGGS, the resolution whose median cell area is closest (in
log-ratio) to the reference H3 res-9 cell — the resolution baked into each
generator's TARGET (and survey.py) so all systems compare cells of roughly the
same size.

Reads the pre-generated *small* cell sets (scripts/dggs/cells/, `just gen-cells`
first) and measures area with `sparea` (spherical-polygon area). No skar, no
DGGS libraries — it runs natively off the Parquet rings. The small set spans
every resolution, which is exactly the scan a calibration needs.

Adding a new DGGS: generate its small set, add a SCAN range below, run this to
get the pick, then bake it into the generator's TARGET + survey.py.

Run with:  just calibrate   (or: uv run scripts/dggs/calibrate.py)
No CLI args (project convention) — edit the constants below in place.
"""

import sys
from pathlib import Path

import numpy as np

import sparea

# _common.py (the cell-set reader) lives in the cells/ subfolder.
sys.path.insert(0, str(Path(__file__).resolve().parent / 'cells'))
import _common as cells  # noqa: E402

# ----- knobs -------------------------------------------------------------
TARGET = ('h3', cells.TARGET_RES['h3'])   # reference system + resolution
# Candidate resolutions to search per system (each within its small set range).
SCAN = {
    's2': range(10, 20),
    'a5': range(8, 20),
    'isea7h': range(0, 16),
    'ivea7h': range(0, 16),
}
# -------------------------------------------------------------------------


def cell_area(dggs, res):
    """Median cell area over the small set's cells, in sparea's native units.

    Only ratios to the reference matter (the pick is by log-ratio), so there's
    no need to convert steradians -> km^2 — the scale factor cancels.
    """
    a = [sparea.area(ring, geo='latlng')
         for _cid, ring in cells.load_cells(dggs, res, 'small')]
    return float(np.median(a))


def main():
    tsys, tres = TARGET
    target = cell_area(tsys, tres)
    print(f'target: {tsys} r{tres} median area = {target:.4e} sr  '
          f'(small set, N={cells.N_SMALL}, seed={cells.SEED:#x})\n')

    for sys, scan in SCAN.items():
        rows = [(res, cell_area(sys, res)) for res in scan]
        best = min(rows, key=lambda r: abs(np.log(r[1] / target)))
        print(f'--- {sys} (target {tsys} r{tres}) ---')
        print(f'{"res":>4} {"area_sr":>11} {"ratio":>8}')
        for res, area in rows:
            mark = '  <== pick' if res == best[0] else ''
            print(f'{res:>4} {area:>11.4e} {area / target:>8.3f}{mark}')
        print(f'-> {sys} r{best[0]}  ({best[1] / target:.3f}x {tsys} r{tres})\n')


if __name__ == '__main__':
    main()

"""Solve aspect ratios for the full-enumeration ivea7h geometry.

Pass two of the full-globe pipeline: reads the binaries written by
gen_ivea7h_full_geom.py (every cell at r5/r6), solves each cell with skar
(native, no DGGS library), and writes one more binary per level:

  ivea7h_r{res}_ar.f32   Float32 aspect ratio per cell (NaN = did-not-converge),
                         len = n_cells. Order matches the _idx.u32 cells.

The browser (globe_full.html) loads pos + idx + ar and colors on the GPU.

Run with:  just web-full   (after `just web-full-geom`)
No CLI args (project convention) — edit RES_LIST below.
"""

import time
from pathlib import Path

import numpy as np

import skar

RES_LIST = [1, 2, 3, 5, 6]
GAP_TOL = 1e-6
OUT = Path(__file__).resolve().parent / 'out' / 'full'


def main():
    for res in RES_LIST:
        pos_path = OUT / f'ivea7h_r{res}_pos.f32'
        idx_path = OUT / f'ivea7h_r{res}_idx.u32'
        ar_path = OUT / f'ivea7h_r{res}_ar.f32'
        if not pos_path.exists():
            print(f'r{res}: {pos_path.name} missing — run `just web-full-geom` first.')
            continue
        if ar_path.exists() and ar_path.stat().st_mtime >= pos_path.stat().st_mtime:
            print(f'r{res}: {ar_path.name} up to date — skipping.')
            continue
        pos = np.fromfile(pos_path, dtype='<f4').reshape(-1, 2)   # [lng, lat]
        idx = np.fromfile(idx_path, dtype='<u4')
        n = len(idx) - 1
        print(f'ivea7h r{res}: solving {n:,} cells...', flush=True)
        t0 = time.perf_counter()
        ar = np.empty(n, dtype='<f4')
        dnc = 0
        for i in range(n):
            ring = pos[idx[i]:idx[i + 1]]                 # (k, 2) [lng, lat]
            latlng = ring[:, ::-1].astype(float)          # (k, 2) [lat, lng]
            r = skar.solve(skar.to_vec3(latlng, geo='latlng_deg'),
                           geo='vec3', gap_tol=GAP_TOL)
            if isinstance(r, skar.Converged):
                ar[i] = r.aspect_ratio
            else:
                ar[i] = np.nan
                dnc += 1
            if (i + 1) % 200_000 == 0:
                print(f'  {i + 1:,}/{n:,}', flush=True)
        ar.tofile(ar_path)
        finite = ar[np.isfinite(ar)]
        print(f'  wrote {ar_path.name}: DNC {dnc:,}, '
              f'AR min {finite.min():.4f} / median {np.median(finite):.4f} / '
              f'max {finite.max():.4f}, {time.perf_counter() - t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()

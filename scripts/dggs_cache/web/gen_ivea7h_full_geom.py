# /// script
# requires-python = "==3.13.*"
# dependencies = ["dggal", "numpy>=1.24"]
# ///
"""Enumerate EVERY cell of SYSTEM at each RES_LIST level and stream its
geometry to binary.

Unlike the sampled Parquet cache (capped at ~100k cells/level), this walks the
full grid — up to 168,072 cells at r5 and 1,176,492 at r6 for ivea7h — for a
complete-coverage globe. DGGAL only (no skar): emits raw geometry; aspect
ratios are solved natively afterwards by build_ivea7h_full_ar.py. Levels whose
binaries already exist are skipped (geometry is deterministic).

Per level it writes two little-endian binaries to web/out/full/ (gitignored):
  {SYSTEM}_r{res}_pos.f32   Float32 [lng, lat] vertex pairs, flattened, open
                            rings as DGGAL yields them — ajglobe converts each
                            vertex to a unit vector, so no antimeridian or
                            winding preprocessing is needed.
  {SYSTEM}_r{res}_idx.u32   Uint32 start indices (len = n_cells + 1), in
                            vertices, marking where each cell's ring begins.
Counts are implicit: n_verts = pos.size/2, n_cells = idx.size - 1.

dggal ships an arch-broken macOS arm64 wheel, so on Apple Silicon this re-execs
under an x86_64 (Rosetta) Python 3.13 — same trick as cells/gen_dggal.py. Run:
    uv run scripts/dggs_cache/web/gen_ivea7h_full_geom.py
No CLI args (project convention) — edit SYSTEM/RES_LIST below.
"""

import ctypes
import glob
import importlib.util
import os
import platform
import sys
import time
from pathlib import Path

SYSTEM = 'ivea7h'          # DGGAL grid to enumerate; its class is SYSTEM.upper()
RES_LIST = [1, 2, 3, 5, 6]

if (sys.platform == 'darwin' and platform.machine() == 'arm64'
        and not os.environ.get('_DGGAL_ROSETTA')):
    os.environ['_DGGAL_ROSETTA'] = '1'
    os.execvp('uv', ['uv', 'run', '--python', 'cpython-3.13-macos-x86_64',
                     os.path.abspath(__file__)])

import numpy as np  # noqa: E402


def _preload_native():
    for pkg, stem in (('ecrt', 'libecrt'), ('dggal', 'libdggal')):
        spec = importlib.util.find_spec(pkg)
        if spec is None or not spec.origin:
            continue
        libdir = os.path.join(os.path.dirname(spec.origin), 'lib')
        for lib in sorted(glob.glob(os.path.join(libdir, stem + '.*'))):
            try:
                ctypes.CDLL(lib, mode=ctypes.RTLD_GLOBAL)
            except OSError:
                pass


try:
    import dggal  # noqa: F401
except ImportError:
    _preload_native()
    import dggal  # noqa: F401

from dggal import *  # noqa: E402,F401,F403

_app = Application(appGlobals=globals())
pydggal_setup(_app)

OUT = Path(__file__).resolve().parent / 'out' / 'full'


def cell_ring(dggrs, zone):
    """Open [(lat, lng), ...] ring for a zone (closing repeat stripped)."""
    ring = [(float(p.lat), float(p.lon)) for p in dggrs.getZoneWGS84Vertices(zone)]
    if len(ring) >= 2 and ring[0] == ring[-1]:
        ring = ring[:-1]
    return ring


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    dggrs = globals()[SYSTEM.upper()]()   # DGGAL class is named after the system
    for res in RES_LIST:
        pos_path = OUT / f'{SYSTEM}_r{res}_pos.f32'
        idx_path = OUT / f'{SYSTEM}_r{res}_idx.u32'
        if pos_path.exists() and idx_path.exists():
            print(f'r{res}: {pos_path.name} exists — skipping (geometry is '
                  f'deterministic; delete to regenerate).', flush=True)
            continue
        n_cells = int(dggrs.countZones(res))
        print(f'{SYSTEM} r{res}: {n_cells:,} cells — enumerating...', flush=True)
        t0 = time.perf_counter()
        starts = [0]
        nverts = 0
        done = 0
        with open(pos_path, 'wb') as fpos:
            for zone in dggrs.listZones(res, wholeWorld):
                ring = cell_ring(dggrs, zone)
                # [lat, lng] rows -> [lng, lat] (the browser axis order)
                arr = np.asarray(ring, dtype='<f4')[:, ::-1]
                fpos.write(arr.tobytes())
                nverts += len(ring)
                starts.append(nverts)
                done += 1
                if done % 100_000 == 0:
                    print(f'  {done:,}/{n_cells:,}', flush=True)
        np.asarray(starts, dtype='<u4').tofile(idx_path)
        mb = (pos_path.stat().st_size + idx_path.stat().st_size) / 1e6
        print(f'  wrote {pos_path.name} + {idx_path.name}: {nverts:,} verts, '
              f'{mb:.1f} MB, {time.perf_counter() - t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()

# /// script
# requires-python = "==3.13.*"
# dependencies = ["dggal", "numpy>=1.24", "pyarrow>=15"]
# ///
"""Generate random DGGAL cell sets (ISEA7H, IVEA7H) -> Parquet.

One file per level (0..the grid's finest) for each system — coarse levels with
<= N cells are enumerated in full, finer ones sample N (schema + write logic in
_common.py). The DGGAL setup + the handful of zone/vertex/count calls cell
generation needs are inlined below; the fuller DGGAL glue for the live analyses
lives in ../../dggs_old/dggal_common.py.

dggal ships an arch-broken macOS arm64 wheel, so on Apple Silicon this script
re-execs itself under an x86_64 (Rosetta) Python 3.13 (where the wheel is
self-consistent) — see the guard below. So it runs the same everywhere:

    uv run scripts/dggs_cache/cells/gen_dggal.py

No CLI args (project convention).
"""

import ctypes
import glob
import importlib.util
import os
import platform
import sys

# On Apple Silicon, dggal's arm64 wheel bundles x86_64 dylibs and won't load, so
# re-exec under an x86_64 (Rosetta) Python 3.13 before importing it — keeping the
# command a plain `uv run` on every platform. Linux wheels are correct (no-op
# there). The env guard breaks the re-exec loop; under Rosetta machine() reports
# 'x86_64', so the condition is false on the second pass anyway.
if (sys.platform == 'darwin' and platform.machine() == 'arm64'
        and not os.environ.get('_DGGAL_ROSETTA')):
    os.environ['_DGGAL_ROSETTA'] = '1'
    os.execvp('uv', ['uv', 'run', '--python', 'cpython-3.13-macos-x86_64',
                     os.path.abspath(__file__)])

import _common  # noqa: E402  (sibling: the cache core)


# ----- DGGAL setup (only what cell generation needs) ---------------------
def _preload_native():
    """dlopen libecrt/libdggal RTLD_GLOBAL so flat-namespace symbols resolve —
    belt-and-suspenders for the arm64 wheel's missing dylib load commands."""
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
    import dggal as _dggal  # noqa: F401
except ImportError:
    _preload_native()
    import dggal as _dggal  # noqa: F401

from dggal import *  # noqa: E402,F401,F403  upstream-documented setup pattern

_app = Application(appGlobals=globals())
pydggal_setup(_app)

# DGGRS class name per system (adding a grid is one row here).
SYSTEMS = {'isea7h': 'ISEA7H', 'ivea7h': 'IVEA7H', 'rhealpix': 'rHEALPix'}


class Adapter:
    """The DGGAL calls cell generation needs, wrapping one DGGRS instance."""

    def __init__(self, cls):
        self.dggrs = globals()[cls]()

    def zone_at(self, level, lng, lat):
        """Zone at `level` containing the (lng, lat) point."""
        return self.dggrs.getZoneFromWGS84Centroid(
            level, GeoPoint(float(lat), float(lng)))

    def cell_boundary(self, zone):
        """Corner vertices as [(lat, lng), ...] deg, open ring (closing repeat
        stripped; handles hexagons and the 12 pentagons)."""
        ring = [(float(p.lat), float(p.lon))
                for p in self.dggrs.getZoneWGS84Vertices(zone)]
        # DGGAL's corner method collapses 2 of 4 corners to (0, 0) for some
        # rHEALPix polar-cap/equatorial-seam cells (a dggal bug). The collapse
        # shows up as duplicate vertices; fall back to the edge-refined boundary,
        # which traces the real cell correctly.
        if len(set(ring)) < len(ring):
            ring = [(float(p.lat), float(p.lon))
                    for p in self.dggrs.getZoneRefinedWGS84Vertices(zone, 0)]
        if len(ring) >= 2 and ring[0] == ring[-1]:
            ring = ring[:-1]
        return ring

    def cid_str(self, zone):
        return self.dggrs.getZoneTextID(zone)

    def count(self, level):
        return int(self.dggrs.countZones(level))

    def max_level(self):
        return self.dggrs.getMaxDGGRSZoneLevel()

    def enumerate(self, level):
        yield from self.dggrs.listZones(level, wholeWorld)


if __name__ == '__main__':
    for name, cls in SYSTEMS.items():
        ad = Adapter(cls)
        _common.generate_levels(
            name, ad.max_level(),
            latlng_to_cell=lambda res, lat, lng, _ad=ad: _ad.zone_at(res, lng, lat),
            cid_str=ad.cid_str,
            cell_boundary=ad.cell_boundary,
            count_at=lambda res, _ad=ad: _ad.count(res),
            enumerate_at=lambda res, _ad=ad: _ad.enumerate(res),
        )

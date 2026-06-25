"""Shared DGGAL binding glue for the DGGS survey scripts.

DGGAL (Ecere's Discrete Global Grid Abstraction Library, `pip install dggal`,
BSD-3-Clause) exposes many DGGRSs through a single `DGGRS` API. This module
initializes the DGGAL `Application` once at import and wraps a DGGRS instance in
the adapter shape survey.py / calibrate.py / dnc_sweep.py expect: count /
enumerate / sample / verts / cid_str (+ area_km2 for calibrate). Each new DGGAL
system is then a one-liner: `Adapter('<DGGRSClass>')` (e.g. 'ISEA7H', 'ISEA3H',
'IVEA7H', 'rHEALPix').

Vertices come back corners-only as an `(M, 3)` unit-vec3 array (M = 6 for
hexagons, 5 for the 12 pentagons), matching the H3/S2/A5 adapters — no repeated
closing vertex.

Platform note: dggal/ecrt only publish an arch-broken macOS arm64 wheel (it
bundles x86_64 dylibs), so the DGGS just-targets run these under an x86_64
(Rosetta) Python 3.13 in a separate env (.venv-dggs) — see the justfile. The
guarded dlopen below is a belt-and-suspenders fallback for that arm64 wheel's
missing dylib load commands; it's a no-op where the wheel is self-consistent.
"""

import ctypes
import glob
import importlib.util
import os

import numpy as np


def _preload_native():
    """dlopen libecrt/libdggal RTLD_GLOBAL so flat-namespace symbols resolve."""
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

# ----- shared constants (mirrors the survey scripts) ---------------------
R_KM = 6371.0088            # mean Earth radius; steradian -> km^2 is R^2
SR2KM2 = R_KM * R_KM
CHUNK = 50_000              # sampling batch size (keeps memory flat)


def sample_uniform_lonlat(n, rng):
    """Uniform-on-sphere samples as (lon_deg, lat_deg), shape (n, 2)."""
    lon = 360.0 * rng.random(n) - 180.0
    lat = np.degrees(np.arcsin(2.0 * rng.random(n) - 1.0))  # equal-area in lat
    return np.column_stack([lon, lat])


def latlng_ring(points):
    """DGGAL WGS84 vertex points (`.lat`/`.lon` Degrees) -> [(lat, lon), ...].

    Corners only: strips a closing repeat if present (matches the H3/S2/A5
    adapters; handles hexagons and the 12 pentagons). Shared by the Adapter
    and validate_corners.py (which feeds it edge-refined vertices instead).
    """
    ring = [(float(p.lat), float(p.lon)) for p in points]
    if len(ring) >= 2 and ring[0] == ring[-1]:
        ring = ring[:-1]
    return ring


class Adapter:
    """Wrap one DGGAL DGGRS in the survey scripts' adapter shape.

    `cls` is the DGGRS class name as exposed by dggal (e.g. 'ISEA7H').
    """

    def __init__(self, cls):
        self.name = cls
        self.dggrs = globals()[cls]()

    # ----- geometry -----------------------------------------------------
    def ring_latlng(self, zone):
        """Corner vertices of `zone` as [(lat, lon), ...] deg, open ring."""
        return latlng_ring(self.dggrs.getZoneWGS84Vertices(zone))

    def verts(self, zone):
        """Corner vertices as an (M, 3) unit-vec3 array (corners only)."""
        import skar  # lazy: the skar-free cell generators import this module too
        return skar.to_vec3(self.ring_latlng(zone), geo='latlng_deg')

    # ----- ids / counts -------------------------------------------------
    def cid_str(self, zone):
        return self.dggrs.getZoneTextID(zone)

    def count(self, level):
        return int(self.dggrs.countZones(level))

    def max_level(self):
        return self.dggrs.getMaxDGGRSZoneLevel()

    # ----- cell streams -------------------------------------------------
    def enumerate(self, level):
        """Every zone at `level`, whole world."""
        yield from self.dggrs.listZones(level, wholeWorld)

    def zone_at(self, level, lon, lat):
        """The zone at `level` containing the (lon, lat) point."""
        return self.dggrs.getZoneFromWGS84Centroid(
            level, GeoPoint(float(lat), float(lon)))

    def sample(self, level, n, rng):
        """`n` zones from uniform-on-sphere points (with repeats)."""
        done = 0
        while done < n:
            k = min(CHUNK, n - done)
            for lon, lat in sample_uniform_lonlat(k, rng):
                yield self.zone_at(level, lon, lat)
            done += k

    def iter_sample(self, level, n, seed):
        """Yield `(cid_str, verts)` for the distinct cells among `n` samples.

        Dedup'd survey stream: draws `n` uniform-on-sphere points, maps each to
        its containing zone, and yields each new zone once. survey.py's
        per-system iterator is a thin wrapper over this.
        """
        rng = np.random.default_rng(seed)
        seen = set()
        for lon, lat in sample_uniform_lonlat(n, rng):
            zone = self.zone_at(level, lon, lat)
            if zone in seen:
                continue
            seen.add(zone)
            yield self.cid_str(zone), self.verts(zone)

    # ----- calibrate ----------------------------------------------------
    def area_km2(self, level, n, seed):
        """Median cell area (km^2) over `n` sampled cells, skar-free."""
        import sparea
        rng = np.random.default_rng(seed)
        a = [sparea.area(self.ring_latlng(self.zone_at(level, lon, lat)),
                         geo='latlng')
             for lon, lat in sample_uniform_lonlat(n, rng)]
        return float(np.median(a)) * SR2KM2


# ----- registered DGGAL grids -------------------------------------------
# One row per DGGAL system in the comparison — the single place to edit to add
# one. `cls` is the DGGRS class name; `res` the H3-r9-matched level (from
# calibrate.py); `scan` the calibrate search range; `color` the matplotlib
# slot. calibrate.py / survey.py / dnc_sweep.py / validate_corners.py each loop
# this dict to register the system, so no per-system functions are needed.
DGGAL_SYSTEMS = {
    'isea7h': dict(cls='ISEA7H', color='C3', res=10, scan=range(0, 16)),
    'ivea7h': dict(cls='IVEA7H', color='C4', res=10, scan=range(0, 16)),
}

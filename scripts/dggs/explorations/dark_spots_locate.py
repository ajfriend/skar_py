"""Locate and characterize the ISEA7H "dark spots" (sharp low-AR cells).

Random sampling at fine resolution rarely lands on the rare low-AR cells, so:
(1) sample a global lon/lat grid, confirm every solve converges (not DNC/nan),
    and list the darkest grid points;
(2) fine-zoom the darkest and walk outward — revealing that these are *razor-
    sharp* dips (AR ~1.0 at a point, ~1.34 just 0.05 deg away), embedded in the
    high-distortion seams (not the smooth face-centroid lows).

See dark_spots_verify.py (not a skar bug) and dark_spots_cells.py (the grid
geometry there). Run under the x86_64 (Rosetta) env — see ../README.md:
    UV_PROJECT_ENVIRONMENT=.venv-dggs uv run --no-sync \
        scripts/dggs/explorations/dark_spots_locate.py
"""

import numpy as np

import skar

from _common import dc

RES = 10
ad = dc.Adapter('ISEA7H')


def ar_at(lon, lat):
    z = ad.zone_at(RES, float(lon), float(lat))
    r = skar.solve(ad.verts(z), geo='vec3')
    ok = isinstance(r, skar.Converged)
    return (r.aspect_ratio if ok else np.nan), ok, (None if ok else r.status)


NLON, NLAT = 720, 360
lons = np.linspace(-180, 180, NLON)
lats = np.linspace(-90, 90, NLAT)
img = np.empty((NLAT, NLON))
nconv = ntot = 0
statuses = {}
for i, lat in enumerate(lats):
    for j, lon in enumerate(lons):
        a, ok, st = ar_at(lon, lat)
        img[i, j] = a
        ntot += 1
        nconv += ok
        if not ok:
            statuses[st] = statuses.get(st, 0) + 1
print(f'converged {nconv}/{ntot}; non-converged {statuses}; '
      f'nan in field {int(np.isnan(img).sum())}')
print(f'field AR range: {np.nanmin(img):.5f} .. {np.nanmax(img):.5f}')

thr = 1.02
ys, xs = np.where(img < thr)
print(f'\n{len(ys)} grid points with AR < {thr}; lowest 25 (lat, lon, AR):')
order = np.argsort(img[ys, xs])
for k in order[:25]:
    y, x = ys[k], xs[k]
    print(f'  lat {lats[y]:6.1f}  lon {lons[x]:7.1f}  AR {img[y, x]:.4f}')

y0, x0 = np.unravel_index(np.nanargmin(img), img.shape)
lat0, lon0 = lats[y0], lons[x0]
b, n = 2.0, 121
zlas = np.linspace(lat0 - b, lat0 + b, n)
zlos = np.linspace(lon0 - b, lon0 + b, n)
zimg = np.empty((n, n))
for i, la in enumerate(zlas):
    for j, lo in enumerate(zlos):
        zimg[i, j], _, _ = ar_at(lo, la)
print(f'\nzoom on darkest (lat {lat0:.1f}, lon {lon0:.1f}), +/-{b} deg:')
print(f'  AR {np.nanmin(zimg):.5f} .. {np.nanmax(zimg):.5f}')
ci = n // 2
print('  AR by ring radius:')
for rdeg in (0.0, 0.25, 0.5, 1.0, 1.5, 2.0):
    di = int(rdeg / b * (n // 2))
    vals = [zimg[ci + dy, ci + dx] for dy, dx in
            ((di, 0), (-di, 0), (0, di), (0, -di))
            if 0 <= ci + dy < n and 0 <= ci + dx < n]
    print(f'    r={rdeg:4.2f} deg: AR ~ {np.nanmean(vals):.4f}')

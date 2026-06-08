"""Spatial map of per-cell AR for ISEA7H vs IVEA7H.

AR ~ the projection's local Tissot anisotropy ratio (resolution-independent for
small cells), so a lon/lat map of cell AR pictures each projection's shape-
distortion field. Two figures:
  ar_heatmap_global.png  — global equirectangular, shared color scale. Shows the
                           20-fold icosahedral pattern: low-distortion face
                           centroids and bright high-distortion seams (ISEA);
                           a near-uniform field (IVEA, no bright seams).
  ar_heatmap_face.png    — gnomonic zoom on one face centroid, per-grid scale.
                           ISEA: low core + sharp 6-pointed seam star; IVEA:
                           smooth bounded radial pinwheel.

Run under the x86_64 (Rosetta) env — see ../README.md "Platform note":
    UV_PROJECT_ENVIRONMENT=.venv-dggs uv run --no-sync \
        scripts/dggs/explorations/ar_heatmap.py
"""

from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from _common import aspect_ratio, dc, tangent_basis

RES = 10
GRIDS = [('ISEA7H', (-58.3971, -168.80)), ('IVEA7H', (58.3971, 11.20))]
OUT = Path(__file__).resolve().parent / 'out'
NG_LON, NG_LAT = 720, 360         # global grid
NF = 500                          # face-zoom grid (NF x NF)
FACE_HALF = 0.7                    # tan(half-angle) ~ 35 deg


def ar_at(ad, lon, lat):
    return aspect_ratio(ad.verts(ad.zone_at(RES, float(lon), float(lat))))


def global_field(ad):
    lons = np.linspace(-180, 180, NG_LON)
    lats = np.linspace(-90, 90, NG_LAT)
    img = np.empty((NG_LAT, NG_LON))
    for i, lat in enumerate(lats):
        for j, lon in enumerate(lons):
            img[i, j] = ar_at(ad, lon, lat)
    return img


def face_field(ad, lat0, lon0):
    c, e1, e2 = tangent_basis(lat0, lon0)
    ts = np.linspace(-FACE_HALF, FACE_HALF, NF)
    img = np.empty((NF, NF))
    for i, y in enumerate(ts):
        for j, x in enumerate(ts):
            d = c + x * e1 + y * e2
            d /= np.linalg.norm(d)
            lat = np.degrees(np.arcsin(d[2]))
            lon = np.degrees(np.arctan2(d[1], d[0]))
            img[i, j] = ar_at(ad, lon, lat)
    return img


OUT.mkdir(parents=True, exist_ok=True)

gfields = {}
for cls, _ in GRIDS:
    print(f'global {cls} ...', flush=True)
    gfields[cls] = global_field(dc.Adapter(cls))
vmax = max(np.nanmax(f) for f in gfields.values())
fig, axes = plt.subplots(len(GRIDS), 1, figsize=(11, 9))
for ax, (cls, _) in zip(axes, GRIDS):
    im = ax.imshow(gfields[cls], origin='lower', extent=[-180, 180, -90, 90],
                   aspect='auto', cmap='viridis', vmin=1.0, vmax=vmax)
    ax.set_title(f'{cls} r{RES} cell aspect ratio (global)', fontsize=10)
    ax.set_ylabel('lat')
    fig.colorbar(im, ax=ax, shrink=0.9, label='AR')
axes[-1].set_xlabel('lon')
fig.suptitle('DGGAL cell aspect-ratio field (shared scale)', fontsize=12)
fig.tight_layout()
fig.savefig(OUT / 'ar_heatmap_global.png', dpi=170)
plt.close(fig)
print('wrote', OUT / 'ar_heatmap_global.png')

ffields = {}
for cls, (lat0, lon0) in GRIDS:
    print(f'face {cls} ...', flush=True)
    ffields[cls] = face_field(dc.Adapter(cls), lat0, lon0)
fig, axes = plt.subplots(1, len(GRIDS), figsize=(12, 5.5))
for ax, (cls, _) in zip(axes, GRIDS):
    f = ffields[cls]
    im = ax.imshow(f, origin='lower', extent=[-FACE_HALF, FACE_HALF] * 2,
                   cmap='viridis')
    ax.set_title(f'{cls} r{RES}  (face centroid; AR {np.nanmin(f):.3f}'
                 f'-{np.nanmax(f):.3f})', fontsize=10)
    ax.set_aspect('equal')
    fig.colorbar(im, ax=ax, shrink=0.85, label='AR')
fig.suptitle('DGGAL cell aspect-ratio over one icosahedral face '
             '(gnomonic, per-grid scale)', fontsize=12)
fig.tight_layout()
fig.savefig(OUT / 'ar_heatmap_face.png', dpi=170)
plt.close(fig)
print('wrote', OUT / 'ar_heatmap_face.png')

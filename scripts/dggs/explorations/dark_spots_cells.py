"""Plot an ISEA7H dark-spot cell + its neighbors (orthographic), with ellipses.

Zoom on a known spike cell. Each cell: refined boundary filled by enclosing-cone
AR, plus its minimum-enclosing ellipse overlaid (black; red for the spike).
Shows that the grid there is a valid continuous tiling of *irregular* hexagons:
the neighbors' bounding ellipses are all elongated along the same seam axis
(AR ~1.34), while the spike cell's is nearly circular (AR ~1.0) — the sharp
shape transition characteristic of a Snyder-ISEA construction seam.

Run under the x86_64 (Rosetta) env — see ../README.md:
    UV_PROJECT_ENVIRONMENT=.venv-dggs uv run --no-sync \
        scripts/dggs/explorations/dark_spots_cells.py
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as MPoly

import skar

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import dggal_common as dc  # noqa: E402  (needs the sys.path insert above)

RES = 10
HALF = 0.007                      # bbox half-size (deg) around the spike
R_KM = 6371.0088
SPIKE_LATLON = (-71.90959, 140.97260)
OUT = Path(__file__).resolve().parent / 'out' / 'dark_spots_cells.png'
ad = dc.Adapter('ISEA7H')


def unit(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo),
                     np.sin(la)])


def refined_xy(z, e1, e2, n=10):
    pts = ad.dggrs.getZoneRefinedWGS84Vertices(z, n)
    ring = [(float(p.lat), float(p.lon)) for p in pts]
    v = skar.to_vec3(ring, geo='latlng_deg')
    return np.column_stack([(v @ e1) * R_KM, (v @ e2) * R_KM])


def corner_xy(z, e1, e2):
    v = ad.verts(z)
    return np.column_stack([(v @ e1) * R_KM, (v @ e2) * R_KM])


def ar(z):
    r = skar.solve(ad.verts(z), geo='vec3')
    return r.aspect_ratio if isinstance(r, skar.Converged) else np.nan


def mvee(P, tol=1e-10):
    N, d = P.shape
    Q = np.vstack([P.T, np.ones(N)])
    u = np.ones(N) / N
    for _ in range(100000):
        X = Q @ np.diag(u) @ Q.T
        M = np.einsum('ij,ji->i', Q.T @ np.linalg.inv(X), Q)
        j = int(np.argmax(M))
        step = (M[j] - d - 1) / ((d + 1) * (M[j] - 1))
        un = (1 - step) * u
        un[j] += step
        if np.linalg.norm(un - u) < tol:
            u = un
            break
        u = un
    c = P.T @ u
    A = np.linalg.inv((P.T @ np.diag(u) @ P) - np.outer(c, c)) / d
    return c, A


def ellipse_pts(c, A, n=120):
    ev, V = np.linalg.eigh(np.linalg.inv(A))
    t = np.linspace(0, 2 * np.pi, n)
    return (V @ np.diag(np.sqrt(ev)) @ np.vstack([np.cos(t), np.sin(t)])).T + c


la0, lo0 = SPIKE_LATLON
spike = ad.zone_at(RES, lo0, la0)
spike_cid = ad.cid_str(spike)
gc = ad.dggrs.getZoneWGS84Centroid(spike)
clat, clon = float(gc.lat), float(gc.lon)
c = unit(clat, clon)
e1 = np.cross([0, 0, 1.0], c)
e1 /= np.linalg.norm(e1)
e2 = np.cross(c, e1)

bbox = dc.GeoExtent((clat - HALF, clon - HALF), (clat + HALF, clon + HALF))
zones = list(ad.dggrs.listZones(RES, bbox))
if not any(ad.cid_str(z) == spike_cid for z in zones):
    zones.append(spike)

polys, ars, ells, spk = [], [], [], []
for z in zones:
    polys.append(refined_xy(z, e1, e2))
    ars.append(ar(z))
    cen, A = mvee(corner_xy(z, e1, e2))
    ells.append(ellipse_pts(cen, A))
    spk.append(ad.cid_str(z) == spike_cid)
ars = np.array(ars)

fig, ax = plt.subplots(figsize=(9.5, 9.5))
pc = PatchCollection([MPoly(p, closed=True) for p in polys],
                     cmap='viridis', alpha=0.55, edgecolor='0.3', linewidths=0.8)
pc.set_array(ars)
pc.set_clim(1.0, np.nanmax(ars))
ax.add_collection(pc)
for p, e, s, a in zip(polys, ells, spk, ars):
    ax.plot(p[:, 0], p[:, 1], '.', color='0.2', ms=2)
    ax.plot(e[:, 0], e[:, 1], '-', color=('red' if s else 'black'),
            lw=(2.4 if s else 1.0))
    ax.text(p[:, 0].mean(), p[:, 1].mean(), f'{a:.2f}', ha='center',
            va='center', fontsize=7, color=('red' if s else '0.1'))
ax.autoscale()
ax.set_aspect('equal')
ax.set_xlabel('east (km)')
ax.set_ylabel('north (km)')
ax.set_title('ISEA7H r10: spike cell (red) + neighbors with min-enclosing '
             'ellipses\n(numbers = enclosing-cone AR)', fontsize=11)
fig.colorbar(pc, ax=ax, shrink=0.85, label='enclosing-cone AR')
fig.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=190)
print(f'{len(zones)} cells; AR {np.nanmin(ars):.3f}..{np.nanmax(ars):.3f}')
print('wrote', OUT)

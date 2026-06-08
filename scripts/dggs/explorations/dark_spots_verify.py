"""Is a dark-spot cell a skar bug, or a real geometry? (independent MVEE check)

skar's enclosing-cone AR should equal the aspect ratio of the minimum-area
enclosing ellipse of the cell's gnomonic-projected vertices. Compute that MVEE
independently (Khachiyan) and compare at a known spike cell. They agree to ~5
decimals (1.00191 vs 1.00191) -> skar is faithful; the spot is a real cell whose
*irregular* hexagon (edge ratio ~1.25, vertex-radius ratio ~1.32) nonetheless
has a near-circular bounding ellipse. Not a bug.

Run under the x86_64 (Rosetta) env — see ../README.md:
    UV_PROJECT_ENVIRONMENT=.venv-dggs uv run --no-sync \
        scripts/dggs/explorations/dark_spots_verify.py
"""

import sys
from pathlib import Path

import numpy as np

import skar

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import dggal_common as dc  # noqa: E402  (needs the sys.path insert above)

RES = 10
ad = dc.Adapter('ISEA7H')


def tangent_xy(z):
    v = ad.verts(z)
    c = v.mean(0)
    c /= np.linalg.norm(c)
    e1 = np.cross([0, 0, 1.0], c)
    if np.linalg.norm(e1) < 1e-9:
        e1 = np.cross([0, 1.0, 0], c)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(c, e1)
    g = v / (v @ c)[:, None]            # gnomonic about centroid
    return np.column_stack([g @ e1, g @ e2])


def mvee_ratio(P, tol=1e-10):
    """Aspect ratio (major/minor) of the min-area enclosing ellipse of P."""
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
    ev = np.linalg.eigvalsh(A)
    return float(np.sqrt(ev.max() / ev.min()))


def skar_ar(z):
    r = skar.solve(ad.verts(z), geo='vec3')
    return r.aspect_ratio if isinstance(r, skar.Converged) else float('nan')


la, lo = -71.90959, 140.97260       # a known ISEA7H spike (from dark_spots_locate)
print(f'{"offset":>7} {"skar":>9} {"MVEE":>9}')
for d in (0.0, 0.05):
    z = ad.zone_at(RES, lo, la + d)
    print(f'{d:7.2f} {skar_ar(z):9.5f} {mvee_ratio(tangent_xy(z)):9.5f}')

z = ad.zone_at(RES, lo, la)
xy = tangent_xy(z)
xy -= xy.mean(0)
edges = np.linalg.norm(np.diff(np.vstack([xy, xy[:1]]), axis=0), axis=1)
radii = np.linalg.norm(xy, axis=1)
print(f'\nspike cell {ad.cid_str(z)}: edge-len ratio {edges.max() / edges.min():.3f}, '
      f'vertex-radius ratio {radii.max() / radii.min():.3f} '
      f'(regular hexagon = 1.000)')
print('=> skar AR matches independent MVEE: not a bug; the irregular hexagon '
      'simply has a near-circular bounding ellipse.')

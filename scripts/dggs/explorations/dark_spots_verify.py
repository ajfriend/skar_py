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

import numpy as np

from _common import (SPIKE_LATLON, aspect_ratio, dc, gnomonic_xy, mvee_ratio,
                     tangent_basis_vec)

RES = 10
ad = dc.Adapter('ISEA7H')


def tangent_xy(z):
    v = ad.verts(z)
    c, e1, e2 = tangent_basis_vec(v.mean(0))    # gnomonic about cell centroid
    return gnomonic_xy(v, c, e1, e2)


la, lo = SPIKE_LATLON
print(f'{"offset":>7} {"skar":>9} {"MVEE":>9}')
for d in (0.0, 0.05):
    z = ad.zone_at(RES, lo, la + d)
    print(f'{d:7.2f} {aspect_ratio(ad.verts(z)):9.5f} {mvee_ratio(tangent_xy(z)):9.5f}')

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

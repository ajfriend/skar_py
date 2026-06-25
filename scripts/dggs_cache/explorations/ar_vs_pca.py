"""Per-cell skar enclosing-cone AR vs PCA second-moment AR (ISEA7H).

Tests whether the near-1.0 tail of the AR distribution is "geometric accidents"
(the minimum-enclosing ellipse collapses to a circle on an *irregular* cell) vs.
the bulk being the real, smooth Tissot / shape-distortion field. Plots a 2-D
histogram of the two metrics per cell and prints the breakdown:
  - bulk cells lie on the skar == PCA diagonal (corr 1.0000, |diff| 0.0000)
    -> the metric IS the real smooth distortion field;
  - the sub-1.0 tail peels off the diagonal (skar ~1.0 while PCA stays
    ~1.1-1.16): an accident of the *grid* -- the ISEA construction makes
    isolated seam cells that are irregular hexagons yet happen to be
    circularly-boundable, which the (correct) enclosing-cone metric reports
    faithfully. Not a metric artifact, and not actually round cells.
(Genuine isotropic cells exist only at the ~20 face centroids and are too rare
by area to be hit by uniform sampling.)

Reads Parquet, no DGGS library, so it runs natively (needs skar built):
    uv run --group cells scripts/dggs_cache/explorations/ar_vs_pca.py
"""

from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

import skar

from _common import aspect_ratio, cells, gnomonic_xy, tangent_basis_vec

RES = cells.TARGET_RES['isea7h']   # = 10, the H3-r9-matched DGGAL level
OUT = Path(__file__).resolve().parent / 'out' / 'ar_vs_pca.png'


def pca_ratio(v):
    """Second-moment (PCA) aspect ratio of a cell's gnomonic-projected verts."""
    c, e1, e2 = tangent_basis_vec(v.mean(0))
    xy = gnomonic_xy(v, c, e1, e2)
    xy -= xy.mean(0)
    s = np.linalg.svd(xy, compute_uv=False)
    return s[0] / s[1]


sk, pc = [], []
for _cid, latlng in cells.load_cells('isea7h', RES):
    v = skar.to_vec3(latlng, geo='latlng_deg')
    sk.append(aspect_ratio(v))
    pc.append(pca_ratio(v))
sk = np.asarray(sk)
pc = np.asarray(pc)

low = sk < 1.05
acc = low & (pc > 1.10)
real = low & (pc < 1.05)
print(f'cells {sk.size:,}')
print(f'skar<1.05: {low.sum()}  -> accidents (PCA>1.10): {acc.sum()}; '
      f'genuinely round (PCA<1.05): {real.sum()}')
mid = (sk > 1.10) & (sk < 1.30)            # bulk
print(f'bulk corr(skar, PCA): {np.corrcoef(sk[mid], pc[mid])[0, 1]:.4f}; '
      f'bulk median |skar-PCA|: {np.median(np.abs(sk[mid] - pc[mid])):.4f}')

fig, ax = plt.subplots(figsize=(8, 8))
h = ax.hist2d(sk, pc, bins=240, norm=mcolors.LogNorm(), cmap='viridis')
lim = [1.0, max(sk.max(), pc.max())]
ax.plot(lim, lim, 'r--', lw=1.2, label='skar = PCA (smooth field)')
ax.set_xlabel('skar enclosing-cone AR (the metric)')
ax.set_ylabel('PCA second-moment AR (interior shape)')
ax.set_title(f'ISEA7H r{RES}: enclosing-cone vs PCA AR, N={sk.size:,}\n'
             'on-diagonal = real smooth distortion; '
             'low-skar/high-PCA = bounding-ellipse accidents')
ax.legend(loc='upper left')
ax.set_aspect('equal')
fig.colorbar(h[3], ax=ax, shrink=0.85, label='cell count (log)')
fig.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=180)
print('wrote', OUT)

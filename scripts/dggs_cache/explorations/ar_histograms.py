"""DGGAL aspect-ratio histograms over the cached cells.

Reads the pre-generated ISEA7H/IVEA7H cell sets at the H3-r9-matched resolution
(scripts/dggs_cache/cells/, `just gen-cells` first) and plots the per-grid AR
distribution (per-grid bins, log y). This is what revealed that ISEA7H is
right-skewed (long tail to ~1.36) while IVEA7H is a flat-topped, sharply-bounded
band (~1.06-1.22).

Reads Parquet, no DGGS library, so it runs natively (needs skar built):
    uv run --group cells scripts/dggs_cache/explorations/ar_histograms.py
"""

from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import skar

from _common import cells

RES = cells.TARGET_RES['isea7h']   # = 10, the H3-r9-matched DGGAL level
BINS = 150
GRIDS = [('isea7h', 'C3'), ('ivea7h', 'C4')]
OUT = Path(__file__).resolve().parent / 'out' / 'ar_histograms.png'

results = {}
for dggs, color in GRIDS:
    ars, dnc = [], 0
    for _cid, latlng in cells.load_cells(dggs, RES):
        r = skar.solve(skar.to_vec3(latlng, geo='latlng_deg'), geo='vec3', gap_tol=1e-6)
        if isinstance(r, skar.Converged):
            ars.append(r.aspect_ratio)
        else:
            dnc += 1
    a = np.asarray(ars)
    results[dggs] = (a, dnc, color)
    deciles = np.percentile(a, np.arange(0, 101, 10))
    print(f'{dggs}: cells={a.size:,} dnc={dnc} min={a.min():.6f} max={a.max():.4f}')
    print('   deciles:', ' '.join(f'{d:.4f}' for d in deciles))

fig, axes = plt.subplots(len(GRIDS), 1, figsize=(9, 8))
for ax, (dggs, (a, dnc, color)) in zip(axes, results.items()):
    bins = np.linspace(1.0, a.max(), BINS + 1)
    ax.hist(a, bins=bins, color=color, edgecolor='white', linewidth=0.2)
    ax.set_yscale('log')
    ax.set_ylabel('count (log)')
    ax.set_xlim(1.0, a.max() * 1.005)
    ax.set_title(f'{dggs.upper()} r{RES}  (N={a.size:,}; min {a.min():.4f}, '
                 f'median {np.median(a):.4f}, max {a.max():.4f}, DNC {dnc})',
                 fontsize=10)
    ax.grid(True, alpha=0.3)
axes[-1].set_xlabel('aspect ratio (per-grid bins, gap_tol = 1e-6)')
fig.suptitle('DGGAL aspect-ratio distributions (~H3 r9 size)', fontsize=12)
fig.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=200)
print('wrote', OUT)

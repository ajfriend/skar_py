"""DGGAL aspect-ratio histograms at a large sample.

Samples N uniform-on-sphere points -> containing cell (dedup'd) for ISEA7H and
IVEA7H at the H3-r9-matched resolution and plots the per-grid AR distribution
(per-grid bins, log y). Larger N => a smoother, more faithful picture of the
area-weighted AR distribution than the survey's N=10k. This is what revealed
that ISEA7H is right-skewed (long tail to ~1.36) while IVEA7H is a flat-topped,
sharply-bounded band (~1.06-1.22). Edit N below (1e6 gives the crispest shape).

Run under the x86_64 (Rosetta) env — see ../README.md "Platform note":
    UV_PROJECT_ENVIRONMENT=.venv-dggs uv run --no-sync \
        scripts/dggs/explorations/ar_histograms.py
"""

from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import skar

from _common import dc

N = 200_000
SEED = 0xC0FFEE
RES = 10
BINS = 150
GRIDS = [('ISEA7H', 'C3'), ('IVEA7H', 'C4')]
OUT = Path(__file__).resolve().parent / 'out' / 'ar_histograms.png'

results = {}
for cls, color in GRIDS:
    ad = dc.Adapter(cls)
    ars, dnc = [], 0
    for _cid, v in ad.iter_sample(RES, N, SEED):
        r = skar.solve(v, geo='vec3', gap_tol=1e-6)
        if isinstance(r, skar.Converged):
            ars.append(r.aspect_ratio)
        else:
            dnc += 1
    a = np.asarray(ars)
    results[cls] = (a, dnc, color)
    deciles = np.percentile(a, np.arange(0, 101, 10))
    print(f'{cls}: cells={a.size:,} dnc={dnc} min={a.min():.6f} max={a.max():.4f}')
    print('   deciles:', ' '.join(f'{d:.4f}' for d in deciles))

fig, axes = plt.subplots(len(GRIDS), 1, figsize=(9, 8))
for ax, (cls, (a, dnc, color)) in zip(axes, results.items()):
    bins = np.linspace(1.0, a.max(), BINS + 1)
    ax.hist(a, bins=bins, color=color, edgecolor='white', linewidth=0.2)
    ax.set_yscale('log')
    ax.set_ylabel('count (log)')
    ax.set_xlim(1.0, a.max() * 1.005)
    ax.set_title(f'{cls} r{RES}  (N={a.size:,}; min {a.min():.4f}, '
                 f'median {np.median(a):.4f}, max {a.max():.4f}, DNC {dnc})',
                 fontsize=10)
    ax.grid(True, alpha=0.3)
axes[-1].set_xlabel('aspect ratio (per-grid bins, gap_tol = 1e-6)')
fig.suptitle(f'DGGAL aspect-ratio distributions, N={N:,} points (~H3 r9 size)',
             fontsize=12)
fig.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=200)
print('wrote', OUT)

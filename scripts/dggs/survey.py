"""DGGS aspect-ratio survey at a commonly-used resolution.

For N random cells at a commonly-used resolution of H3, S2, A5, and the
DGGAL ISEA7H/IVEA7H hex grids, compute the tightest enclosing-cone aspect
ratio with `skar`, then plot the per-system distribution and the
best/worst cell with its enclosing ellipse.

Resolution choice: H3 res 9 is the reference (~0.1 km^2, ~174 m edge —
a typical working resolution, not the metre-scale finest). The others are
set to the resolution whose cell area is closest to an H3 r9 cell:
S2 level 15 (0.76x H3 r9 area), A5 resolution 14 (1.15x), ISEA7H/IVEA7H
r10 (1.65x). They refine in coarser steps than x1, so none lands exactly
on target; these are the nearest in log-area. Recompute with
`just calibrate` (calibrate.py) when adding a new DGGS.

Cells come from Parquet, not sampled here: scripts/dggs/cells/gen_<dggs>.py
(standalone uv-run scripts, each with its own DGGS-library dependency)
write distinct random cells to scripts/dggs/cells/out/. This survey reads
those rings back and solves each — so it imports no DGGS library and runs
natively (no Rosetta). Generate the cell sets first with `just gen-cells`.
Each cell is solved as it streams in; only running aggregates are kept —
the aspect ratios (one float per cell) and the two extreme cells per
system. Nothing new is materialized to disk except the final PNGs.

Solve tolerance: skar's strict gap_tol = 1e-6 default. Every cell at
these resolutions converges at 1e-6. (A band of H3 resolutions, r7-r10,
used to stall at ~1.7e-6 and needed a relaxed 1e-5; skar_zig v0.2.0 fixed
it by lowering the certificate active-set cutoff ACTIVE_THRESH from 1e-6
to 1e-12 — see h3_gap_floor_report.md at repo root.)

Run with:  just dggs   (reads the Parquet cell sets; `just gen-cells` first)
No CLI args (project convention) — edit the constants below in place.
"""

import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import skar

# Cell sets are pre-generated to Parquet by scripts/dggs/cells/gen_<dggs>.py —
# each a standalone uv-run script carrying its own DGGS-library dependency. This
# survey is itself DGGS-library-free: it reads the rings back and solves them.
# Generate them first with `just gen-cells` (see scripts/dggs/cells/README.md).
sys.path.insert(0, str(Path(__file__).resolve().parent / 'cells'))
import _common as cells  # noqa: E402

# ----- knobs -------------------------------------------------------------
# N and SEED must match what the generators wrote (their file-naming constants).
N = 100_000
SEED = 0xC0FFEE
GAP_TOL = 1e-6

# Resolution per system, matched to H3 r9 cell area (median over random cells,
# via calibrate.py): H3 r9 ~0.110 km^2 (target); S2 L15 0.083 km^2 (0.76x);
# A5 r14 0.127 km^2 (1.15x); ISEA7H/IVEA7H r10 0.181 km^2 (1.65x — aperture-7
# steps by 7x, so r10 is the nearest level). Recompute with `just calibrate`.
RES = {'h3': 9, 's2': 15, 'a5': 14, 'isea7h': 10, 'ivea7h': 10}

SYSTEMS = list(RES)
SYS_LABEL = {'h3': 'H3 r9', 's2': 'S2 L15', 'a5': 'A5 r14',
             'isea7h': 'ISEA7H r10', 'ivea7h': 'IVEA7H r10'}
SYS_COLOR = {'h3': 'C0', 's2': 'C1', 'a5': 'C2', 'isea7h': 'C3', 'ivea7h': 'C4'}

OUT_DIR = Path(__file__).resolve().parent / 'out'
N_BINS = 60
DPI = 200
# -------------------------------------------------------------------------


# ----- per-system cell stream: yield (id, (M, 3) unit-vertex array) -------
def iter_cells(name):
    """Stream (id, (M, 3) unit-vertex array) for `name` from its Parquet set."""
    for cid, latlng in cells.load_cells(name, RES[name], N, SEED):
        yield cid, skar.to_vec3(latlng, geo='latlng_deg')


def run_system(name):
    """Stream every cell, solve it, keep aspect ratios + the two extremes."""
    ars = []
    dnc = 0
    best = worst = None
    for cid, verts in iter_cells(name):
        r = skar.solve(verts, geo='vec3', gap_tol=GAP_TOL)
        if not isinstance(r, skar.Converged):
            dnc += 1
            continue
        ar = r.aspect_ratio
        ars.append(ar)
        rec = {'ar': ar, 'id': cid, 'verts': verts, 'result': r}
        if best is None or ar < best['ar']:
            best = rec
        if worst is None or ar > worst['ar']:
            worst = rec
    return {'ars': np.asarray(ars), 'dnc': dnc, 'best': best, 'worst': worst}


# ----- plotting ----------------------------------------------------------
def plot_histograms(results):
    ars = {s: results[s]['ars'] for s in SYSTEMS}
    dnc = {s: results[s]['dnc'] for s in SYSTEMS}

    def stat(a):
        return dict(n=a.size, min=float(a.min()), median=float(np.median(a)),
                    p99=float(np.percentile(a, 99)), max=float(a.max()))

    st = {s: stat(ars[s]) for s in SYSTEMS}

    print(f'{"sys":5} {"n_conv":>8} {"n_dnc":>7} {"min":>10} {"median":>10} {"p99":>10} {"max":>10}')
    for s in SYSTEMS:
        d = st[s]
        print(f'{s:5} {d["n"]:>8} {dnc[s]:>7} {d["min"]:>10.6f} '
              f'{d["median"]:>10.6f} {d["p99"]:>10.6f} {d["max"]:>10.6f}')

    # Shared bins across systems (AR >= 1 by definition) for comparability.
    bins = np.linspace(1.0, max(a.max() for a in ars.values()), N_BINS + 1)

    fig, axes = plt.subplots(len(SYSTEMS), 1, figsize=(8, 9), sharex=True)
    for ax, s in zip(axes, SYSTEMS):
        d = st[s]
        ax.hist(ars[s], bins=bins, color=SYS_COLOR[s], edgecolor='white', linewidth=0.3)
        ax.set_yscale('log')  # AR clusters near the low end with a thin tail
        ax.set_ylabel('count (log)')
        ax.set_title(f'{SYS_LABEL[s]}  (median {d["median"]:.4f}, max {d["max"]:.4f}, DNC {dnc[s]})',
                     fontsize=10)
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel('aspect ratio (shared bins, gap_tol = 1e-6)')
    fig.suptitle('DGGS aspect-ratio distributions (~H3 r9 cell size)', fontsize=12)
    fig.tight_layout()
    out = OUT_DIR / 'histograms.png'
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f'wrote {out}')


def draw_cell(ax, rec, color):
    """Draw a cell's boundary + enclosing ellipse, major axis horizontal."""
    # up=None: no north-up flip (these cells aren't oriented to geographic north).
    xy, semi = skar.project_to_cone(rec['result'], rec['verts'], up=None)
    ring = np.vstack([xy, xy[:1]])
    t = np.linspace(0.0, 2.0 * np.pi, 400)
    ax.plot(ring[:, 0], ring[:, 1], '-o', color=color, lw=1.3, ms=4, label='cell')
    ax.plot(semi[0] * np.cos(t), semi[1] * np.sin(t), '-', color='0.25', lw=1.5,
            label='enclosing ellipse')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.text(0.03, 0.95, f'AR {rec["ar"]:.4f}\nid {rec["id"]}', transform=ax.transAxes,
            va='top', ha='left', fontsize=8,
            bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.85))


def plot_extremes(results):
    fig, axes = plt.subplots(len(SYSTEMS), 2, figsize=(11, 4.3 * len(SYSTEMS)))
    axes[0, 0].set_title('best AR (most circular)', fontsize=12, pad=10)
    axes[0, 1].set_title('worst AR', fontsize=12, pad=10)
    for row, s in enumerate(SYSTEMS):
        for col, kind in ((0, 'best'), (1, 'worst')):
            ax = axes[row, col]
            draw_cell(ax, results[s][kind], SYS_COLOR[s])
            ax.set_xlabel('major axis (m)')
            if col == 0:
                ax.set_ylabel(f'{SYS_LABEL[s]}\nminor axis (m)')
    axes[0, 0].legend(loc='lower right', fontsize=8)
    fig.suptitle('DGGS cells (~H3 r9 cell size): best vs worst aspect ratio\n'
                 '(enclosing-cone cross-section ‖Ax‖ <= b·x; major axis horizontal)',
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = OUT_DIR / 'extremes.png'
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f'wrote {out}')
    for s in SYSTEMS:
        print(f'  {s}: best AR {results[s]["best"]["ar"]:.4f}  ·  '
              f'worst AR {results[s]["worst"]["ar"]:.4f}')


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {}
    for s in SYSTEMS:
        t0 = time.perf_counter()
        results[s] = run_system(s)
        dt = time.perf_counter() - t0
        n = results[s]['ars'].size
        print(f'[{s}] {n} cells solved in {dt:.2f}s (DNC {results[s]["dnc"]})')

    plot_histograms(results)
    plot_extremes(results)


if __name__ == '__main__':
    main()

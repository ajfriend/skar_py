"""DGGS finest-resolution aspect-ratio survey.

For N random cells at the finest resolution of H3, S2, and A5, compute
the tightest enclosing-cone aspect ratio with `skar`, then plot the
per-system distribution and the best/worst cell with its enclosing
ellipse.

Single-pass and file-free: a generator streams one cell at a time
(`(id, unit-vertex array)`), each is solved immediately, and only the
running aggregates are kept — the aspect ratios (one float per cell) and
the two extreme cells per system. Nothing is materialized to disk except
the final PNGs. This is the Python port of the old gen -> JSON -> Zig ->
JSON -> plot pipeline, collapsed now that `skar.solve` is callable here.

Solve tolerance: gap_tol = 1e-3, not skar's strict 1e-6 default. At
finest resolution the S2/A5 cells are sub-metre scatters at an O(1)
point on the sphere (kappa(A) ~ 1e9), so the duality gap floors at
~1e-4-1e-3 and many cells would return `did_not_converge` at 1e-6 — yet
their aspect ratios are accurate regardless (input-precision-limited).
Solving at 1e-3 lets every cell converge so the distribution is complete.

Run with:  just dggs        (or: uv run --group dggs scripts/dggs/survey.py)
No CLI args (project convention) — edit the constants below in place.
"""

import time
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import a5
import h3
import s2sphere

import skar

# ----- knobs -------------------------------------------------------------
N = 10_000
SEED = 0xC0FFEE
GAP_TOL = 1e-3

H3_RES = 15                 # h3 supports 0..15
S2_LEVEL = 30               # s2sphere supports 0..30
A5_RES = a5.MAX_RESOLUTION  # currently 30

OUT_DIR = Path(__file__).resolve().parent / 'out'
N_BINS = 60
DPI = 200

SYSTEMS = ['h3', 's2', 'a5']
SYS_LABEL = {'h3': 'H3 r15', 's2': 'S2 L30', 'a5': 'A5 r30'}
SYS_COLOR = {'h3': 'C0', 's2': 'C1', 'a5': 'C2'}
# -------------------------------------------------------------------------


def sample_uniform_lonlat(n, rng):
    """Uniform-on-sphere samples as (lon_deg, lat_deg), shape (n, 2)."""
    lon = 360.0 * rng.random(n) - 180.0
    lat = np.degrees(np.arcsin(2.0 * rng.random(n) - 1.0))  # equal-area in lat
    return np.column_stack([lon, lat])


# ----- per-system cell streams: yield (id, (M, 3) unit-vertex array) ------
def iter_h3(n, seed):
    rng = np.random.default_rng(seed)
    seen = set()
    for lon, lat in sample_uniform_lonlat(n, rng):
        cid = h3.latlng_to_cell(float(lat), float(lon), H3_RES)
        if cid in seen:
            continue
        seen.add(cid)
        boundary = h3.cell_to_boundary(cid)  # [(lat, lng), ...]
        yield cid, skar.to_vec3(boundary, geo='latlng')


def iter_s2(n, seed):
    rng = np.random.default_rng(seed)
    seen = set()
    for lon, lat in sample_uniform_lonlat(n, rng):
        cid = s2sphere.CellId.from_lat_lng(s2sphere.LatLng.from_degrees(float(lat), float(lon)))
        if S2_LEVEL != 30:
            cid = cid.parent(S2_LEVEL)
        key = cid.id()
        if key in seen:
            continue
        seen.add(key)
        cell = s2sphere.Cell(cid)
        # get_vertex returns a Point that is not necessarily unit length.
        v = np.array([tuple(cell.get_vertex(i))[:3] for i in range(4)], dtype=float)
        v /= np.linalg.norm(v, axis=1, keepdims=True)
        yield format(key, '016x'), v


def iter_a5(n, seed):
    rng = np.random.default_rng(seed)
    seen = set()
    for lon, lat in sample_uniform_lonlat(n, rng):
        cid = a5.lonlat_to_cell((float(lon), float(lat)), A5_RES)
        if cid in seen:
            continue
        seen.add(cid)
        ring = a5.cell_to_boundary(cid)  # closed ring of (lon, lat)
        if len(ring) >= 2 and tuple(ring[0]) == tuple(ring[-1]):
            ring = ring[:-1]
        latlng = [(lat_, lon_) for lon_, lat_ in ring]
        yield a5.u64_to_hex(cid), skar.to_vec3(latlng, geo='latlng_deg')


ITERATORS = {'h3': iter_h3, 's2': iter_s2, 'a5': iter_a5}


def run_system(name):
    """Stream every cell, solve it, keep aspect ratios + the two extremes."""
    ars = []
    dnc = 0
    best = worst = None
    for cid, verts in ITERATORS[name](N, SEED):
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
    axes[-1].set_xlabel('aspect ratio (shared bins, gap_tol = 1e-3)')
    fig.suptitle('DGGS finest-resolution aspect-ratio distributions', fontsize=12)
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
    fig, axes = plt.subplots(3, 2, figsize=(11, 13))
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
    fig.suptitle('DGGS finest-resolution cells: best vs worst aspect ratio\n'
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

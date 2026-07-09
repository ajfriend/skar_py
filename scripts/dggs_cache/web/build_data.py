"""Build the static data for the DGGS aspect-ratio web viewer.

A third reader of the Parquet cell cache (alongside survey.py / dnc_check.py):
reads the pre-generated rings (`just gen-cells` first), solves every cell with
`skar`, and emits browser-friendly JSON + flat binaries into web/out/
(gitignored). No DGGS library, native arch — same decoupling as the rest of
dggs_cache.

Two products:

- out/histograms.json — for every (system, resolution): a FIXED fine-bin
  histogram (NBINS bins over [1, AMAX] + one overflow bin) plus summary stats.
  The fixed grid is what lets the page re-aggregate to any coarser bin width in
  the browser; AMAX is a high global percentile so the working resolutions are
  never pushed into overflow. The page picks the *displayed* x-domain
  dynamically from whatever series are selected.

- out/globe/{sys}_r{res}_{pos.f32,idx.u32,ar.f32,ids.json} — for the coarse
  resolutions only (largest res per system with <= GLOBE_MAX_CELLS cells, and
  everything below it): ajglobe's native flat-binary polygon format. pos is
  Float32 [lng, lat] vertex pairs (open rings, any winding — ajglobe
  triangulates by topology), idx is Uint32 ring start indices
  (len = n_cells + 1), ar is Float32 per cell (NaN = did-not-converge), and
  ids.json is the cell-id strings for hover. Drives the orthographic globes.

- out/manifest.json — what exists (hist + globe resolutions per system), the
  per-system web colors, the fixed bin grid, and the shared globe AR range.

Run with:  just web-data   (`just gen-cells` first)
No CLI args (project convention) — edit the constants below in place.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from matplotlib.colors import to_hex

import skar

# Cell sets are pre-generated to Parquet by scripts/dggs_cache/cells/gen_<dggs>.py.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'cells'))
import _common as cells  # noqa: E402

# ----- knobs -------------------------------------------------------------
GAP_TOL = 1e-6           # solve tolerance (matches survey.py)
NBINS = 256              # fixed fine bins for the stored histograms
AMAX_PCT = 99.99         # global percentile that sets the fixed grid's top edge
GLOBE_MAX_CELLS = 80_000  # a system's globe shows the largest res at/under this

SYSTEMS = list(cells.TARGET_RES)
# matplotlib resolves 'C0'..'C5' itself, so the web colors match the survey
# PNGs exactly.
SYS_COLOR = {s: to_hex(cells.SYS_COLOR[s]) for s in SYSTEMS}
SYS_LABEL = {s: s.upper() for s in SYSTEMS}
RES_PREFIX = {s: 'L' if s == 's2' else 'r' for s in SYSTEMS}  # S2 calls them levels

OUT_DIR = Path(__file__).resolve().parent / 'out'
GLOBE_DIR = OUT_DIR / 'globe'
# -------------------------------------------------------------------------


def solve_ar(latlng):
    """AR for one cell's (M,2) lat/lng ring, or NaN if it didn't converge."""
    r = skar.solve(skar.to_vec3(latlng, geo='latlng_deg'), geo='vec3', gap_tol=GAP_TOL)
    return r.aspect_ratio if isinstance(r, skar.Converged) else np.nan


def sweep_all():
    """Solve every cell of every cached (system, resolution). Returns
    {system: {res: np.ndarray of per-cell ARs}} in load_cells order, NaN where
    the cell didn't converge — build_globe reuses these arrays cell-for-cell,
    so nothing is solved twice."""
    ars = {}
    for s in SYSTEMS:
        ars[s] = {}
        t0 = time.perf_counter()
        total = 0
        for res in cells.available_resolutions(s):
            ars[s][res] = np.asarray([solve_ar(latlng)
                                      for _cid, latlng in cells.load_cells(s, res)])
            total += ars[s][res].size
        print(f'[{s}] {total:,} cells over {len(ars[s])} resolutions '
              f'in {time.perf_counter() - t0:.1f}s')
    return ars


def build_histograms(ars):
    """Fixed-grid histograms + stats per (system, res). The grid is shared
    across every system so the page can overlay and re-bin them on one axis."""
    allars = np.concatenate([a[~np.isnan(a)] for byres in ars.values()
                             for a in byres.values() if a.size])
    amax = float(np.percentile(allars, AMAX_PCT))
    edges = np.linspace(1.0, amax, NBINS + 1)

    data = {}
    for s in SYSTEMS:
        data[s] = {}
        for res, cell_ars in ars[s].items():
            a = cell_ars[~np.isnan(cell_ars)]
            if not a.size:
                continue
            counts = np.histogram(a, bins=edges)[0]
            data[s][str(res)] = {
                'counts': counts.astype(int).tolist(),
                'n': int(a.size),
                'dnc': int(np.isnan(cell_ars).sum()),
                'min': float(a.min()),
                'median': float(np.median(a)),
                'p99': float(np.percentile(a, 99)),
                'max': float(a.max()),
            }
    return {
        'edges': [float(x) for x in edges],   # NBINS+1 edges; counts[i] in [edges[i], edges[i+1])
        'nbins': NBINS,
        'amax': amax,
        'data': data,
    }


def globe_resolutions(s):
    """The coarse resolutions to render for system `s`: the contiguous prefix
    from the base grid up to the last res with <= GLOBE_MAX_CELLS cells. Stops at
    the FIRST over-cap res rather than filtering all of them — past the system's
    target resolution the cache only keeps N_SMALL *sampled* cells, so those deep
    levels also fall under the cap but are sparse scatter, not globe coverage."""
    out = []
    for res in cells.available_resolutions(s):
        nrows = pq.ParquetFile(cells.cells_path(s, res)).metadata.num_rows
        if nrows > GLOBE_MAX_CELLS:
            break
        out.append(res)
    return out


def build_globe(ars):
    """Write ajglobe's flat binaries per coarse (system, res): pos.f32
    ([lng, lat] vertex pairs, open rings), idx.u32 (ring starts), ar.f32
    (NaN = DNC), ids.json. Reuses the swept AR arrays (same load_cells order),
    so this pass only re-reads geometry. Returns {system: [res, ...]} and the
    shared AR max over all globe cells."""
    GLOBE_DIR.mkdir(parents=True, exist_ok=True)
    avail, globe_max = {}, 1.0
    for s in SYSTEMS:
        res_list = globe_resolutions(s)
        avail[s] = res_list
        for res in res_list:
            cell_ars = ars[s][res]
            globe_max = max(globe_max, float(np.nanmax(cell_ars)))
            pos, starts, ids = [], [0], []
            for cid, latlng in cells.load_cells(s, res):
                pos.append(np.asarray(latlng, dtype='<f4')[:, ::-1])  # -> [lng, lat]
                starts.append(starts[-1] + len(latlng))
                ids.append(cid)
            stem = GLOBE_DIR / f'{s}_r{res}'
            np.concatenate(pos).tofile(f'{stem}_pos.f32')
            np.asarray(starts, dtype='<u4').tofile(f'{stem}_idx.u32')
            cell_ars.astype('<f4').tofile(f'{stem}_ar.f32')
            Path(f'{stem}_ids.json').write_text(json.dumps(ids))
            print(f'  globe {s} r{res}: {len(ids)} cells -> {stem.name}_*')
    return avail, globe_max


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print('solving all cells...')
    ars = sweep_all()

    print('building histograms...')
    hist = build_histograms(ars)
    (OUT_DIR / 'histograms.json').write_text(json.dumps(hist))
    print(f'wrote histograms.json (amax={hist["amax"]:.4f}, {NBINS} bins)')

    print('building globe binaries...')
    globe_avail, globe_max = build_globe(ars)

    manifest = {
        'systems': SYSTEMS,
        'colors': SYS_COLOR,
        'labels': SYS_LABEL,
        'res_prefix': RES_PREFIX,
        'target_res': cells.TARGET_RES,
        'hist_res': {s: sorted(int(r) for r in hist['data'].get(s, {})) for s in SYSTEMS},
        'globe_res': globe_avail,
        'globe_ar_max': globe_max,
        'gap_tol': GAP_TOL,
    }
    (OUT_DIR / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(f'wrote manifest.json (globe AR range 1..{globe_max:.4f})')
    print('done.')


if __name__ == '__main__':
    main()

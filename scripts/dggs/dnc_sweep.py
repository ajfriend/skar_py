"""H3/S2/A5 cross-resolution DNC sweep — boundary + non-monotonicity.

Solves cells across every resolution of H3 (0..15), S2 (levels 0..30) and
A5 (resolutions 0..30) with skar's DEFAULT settings — `skar.solve(v,
geo='vec3')`, gap_tol=1e-6, max_outer=100 — and characterizes where
`did_not_converge` (DNC) appears.

H3 stays clean across all its resolutions (its finest, r15, is still ~0.9
m^2). S2/A5 *will* DNC at their finest resolutions: their sub-metre cells
hit a genuine f64 duality-gap floor (kappa ~ 1e9), which is correct
behaviour (see h3_gap_floor_report.md / skar_zig dggs_dnc_test.zig: ~22%
of S2 L30 and ~47% of A5 r30 DNC at 1e-6). What we're hunting for is:
  1. the DNC boundary (onset resolution and how the fraction rises);
  2. NON-MONOTONIC behaviour — a coarser resolution with more DNC than a
     finer one, or a "DNC island" (DNC>0 with a finer resolution at 0),
     i.e. the H3-r7-r10-style surprise that v0.2.0 fixed;
  3. any unexpected DNCs, dumped reproducibly to out/dnc_sweep_cells.txt.

Per resolution: enumerate every cell when there are few enough
(<= ENUMERATE_MAX, exact/noise-free), else sample N_PER_RES random cells
(H3 always also tests all 12 pentagons, even when sampling).

Writes a 2-panel PNG (DNC fraction + worst converged gap vs resolution)
to out/dnc_sweep.png and prints per-system tables + a monotonicity report.

Run with:  just dnc-sweep
       (or: uv run --group dggs scripts/dggs/dnc_sweep.py)
No CLI args (project convention) — edit the constants below in place.
"""

import time
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import a5_fast as a5  # Rust/PyO3 A5 binding (~30x faster than pure-Python pya5)
import h3
import s2sphere

import skar

import dggal_common  # Ecere DGGAL binding glue (ISEA/IVEA hex, rHEALPix, ...)

# ----- knobs -------------------------------------------------------------
SEED = 0xC0FFEE
# random cells per sampled resolution
N_PER_RES = {'h3': 500_000, 's2': 500_000, 'a5': 100_000, 'isea7h': 100_000}
ENUMERATE_MAX = 400_000   # resolutions with <= this many cells: test all
CHUNK = 50_000            # sampling batch size (keeps memory flat)
NOISE_TOL = 5e-4          # monotonicity: ignore DNC-fraction dips below this
MAX_DUMP_PER_RES = 50     # cells written per flagged resolution

OUT_DIR = Path(__file__).resolve().parent / 'out'
PNG = OUT_DIR / 'dnc_sweep.png'
CELLS_FILE = OUT_DIR / 'dnc_sweep_cells.txt'

SYS_COLOR = {'h3': 'C0', 's2': 'C1', 'a5': 'C2', 'isea7h': 'C3'}
SYS_LABEL = {'h3': 'H3', 's2': 'S2', 'a5': 'A5', 'isea7h': 'ISEA7H'}
# -------------------------------------------------------------------------


def sample_uniform_lonlat(n, rng):
    """Uniform-on-sphere samples as (lon_deg, lat_deg), shape (n, 2)."""
    lon = 360.0 * rng.random(n) - 180.0
    lat = np.degrees(np.arcsin(2.0 * rng.random(n) - 1.0))  # equal-area in lat
    return np.column_stack([lon, lat])


# ----- H3 adapter --------------------------------------------------------
def h3_count(res):
    return h3.get_num_cells(res)


def h3_enumerate(res):
    if res == 0:
        yield from h3.get_res0_cells()
    else:
        for c0 in h3.get_res0_cells():
            yield from h3.cell_to_children(c0, res)


def h3_sample(res, n, rng):
    yield from h3.get_pentagons(res)  # always test pentagons, even when sampling
    done = 0
    while done < n:
        k = min(CHUNK, n - done)
        for lon, lat in sample_uniform_lonlat(k, rng):
            yield h3.latlng_to_cell(float(lat), float(lon), res)
        done += k


def h3_verts(cid):
    return skar.to_vec3(h3.cell_to_boundary(cid), geo='latlng')  # ring is (lat, lng)


def h3_id(cid):
    return str(cid)


# ----- S2 adapter --------------------------------------------------------
def s2_count(res):
    return 6 * 4 ** res


def s2_enumerate(res):
    return s2sphere.CellId.walk(res)


def s2_sample(res, n, rng):
    done = 0
    while done < n:
        k = min(CHUNK, n - done)
        for lon, lat in sample_uniform_lonlat(k, rng):
            cid = s2sphere.CellId.from_lat_lng(
                s2sphere.LatLng.from_degrees(float(lat), float(lon)))
            yield cid.parent(res)
        done += k


def s2_verts(cid):
    cell = s2sphere.Cell(cid)
    v = np.array([tuple(cell.get_vertex(i))[:3] for i in range(4)], dtype=float)
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    return v


def s2_id(cid):
    return format(cid.id(), '016x')


# ----- A5 adapter --------------------------------------------------------
def a5_count(res):
    return a5.get_num_cells(res)


def a5_enumerate(res):
    res0 = a5.get_res0_cells()
    if res == 0:
        yield from res0
    else:
        for c0 in res0:
            yield from a5.cell_to_children(c0, res)


def a5_sample(res, n, rng):
    done = 0
    while done < n:
        k = min(CHUNK, n - done)
        for lon, lat in sample_uniform_lonlat(k, rng):
            yield a5.lonlat_to_cell(float(lon), float(lat), res)
        done += k


def a5_verts(cid):
    ring = a5.cell_to_boundary(cid)  # closed ring of (lon, lat)
    if len(ring) >= 2 and tuple(ring[0]) == tuple(ring[-1]):
        ring = ring[:-1]
    latlng = [(lat_, lon_) for lon_, lat_ in ring]
    return skar.to_vec3(latlng, geo='latlng_deg')


def a5_id(cid):
    return a5.u64_to_hex(cid)


# ----- ISEA7H adapter (DGGAL) -------------------------------------------
# The Adapter already exposes count/enumerate/sample/verts/cid_str in the
# right shape, so the SYSTEMS entry references its bound methods directly —
# no per-system wrapper functions (each future DGGAL grid is one dict entry).
_isea7h = dggal_common.Adapter('ISEA7H')


SYSTEMS = {
    'h3': dict(count=h3_count, enumerate=h3_enumerate, sample=h3_sample,
               verts=h3_verts, cid_str=h3_id, res_range=range(0, 16)),
    's2': dict(count=s2_count, enumerate=s2_enumerate, sample=s2_sample,
               verts=s2_verts, cid_str=s2_id, res_range=range(0, 31)),
    'a5': dict(count=a5_count, enumerate=a5_enumerate, sample=a5_sample,
               verts=a5_verts, cid_str=a5_id, res_range=range(0, 31)),
    'isea7h': dict(count=_isea7h.count, enumerate=_isea7h.enumerate,
                   sample=_isea7h.sample, verts=_isea7h.verts,
                   cid_str=_isea7h.cid_str,
                   res_range=range(0, 20)),  # isea7h max level 19
}


# ----- sweep -------------------------------------------------------------
def cells_for_res(sys, res, rng):
    if SYSTEMS[sys]['count'](res) <= ENUMERATE_MAX:
        return 'all', SYSTEMS[sys]['enumerate'](res)
    return 'sample', SYSTEMS[sys]['sample'](res, N_PER_RES[sys], rng)


def sweep_system(sys):
    """Return a list of per-resolution record dicts."""
    adapter = SYSTEMS[sys]
    rng = np.random.default_rng(SEED)
    rows = []
    print(f'\n=== {SYS_LABEL[sys]} ===')
    print(f'{"res":>3} {"mode":>7} {"tested":>10} {"dnc":>8} {"dnc%":>7} '
          f'{"conv_gap":>10} {"dnc_gap":>10} {"max_it":>6} {"secs":>6}')
    for res in adapter['res_range']:
        mode, stream = cells_for_res(sys, res, rng)
        t0 = time.perf_counter()
        tested = dnc = infeas = raised = 0
        conv_worst = 0.0
        max_it = 0
        dnc_gaps = []
        dump = []  # (cid_str, gap, verts) for reproduction
        for cid in stream:
            tested += 1
            try:
                v = adapter['verts'](cid)
                r = skar.solve(v, geo='vec3')  # DEFAULT settings
            except Exception as e:
                raised += 1
                if len(dump) < MAX_DUMP_PER_RES:
                    dump.append((adapter['cid_str'](cid),
                                 f'raised:{type(e).__name__}', None))
                continue
            if isinstance(r, skar.Converged):
                conv_worst = max(conv_worst, r.gap)
                max_it = max(max_it, r.outer_iters)
            elif r.status == 'infeasible':
                infeas += 1
            else:  # did_not_converge
                dnc += 1
                dnc_gaps.append(r.gap)
                if len(dump) < MAX_DUMP_PER_RES:
                    try:
                        verts = adapter['verts'](cid)
                    except Exception:
                        verts = None
                    dump.append((adapter['cid_str'](cid), r.gap, verts))
        dt = time.perf_counter() - t0
        bad = dnc + infeas + raised
        frac = bad / tested if tested else 0.0
        dnc_gap_med = float(np.median(dnc_gaps)) if dnc_gaps else float('nan')
        rows.append(dict(sys=sys, res=res, mode=mode, tested=tested, dnc=dnc,
                         infeas=infeas, raised=raised, bad=bad, frac=frac,
                         conv_worst=conv_worst, dnc_gap_med=dnc_gap_med,
                         max_it=max_it, dump=dump))
        print(f'{res:>3} {mode:>7} {tested:>10,} {bad:>8} {100*frac:>6.2f}% '
              f'{conv_worst:>10.2e} {dnc_gap_med:>10.2e} {max_it:>6} {dt:>6.1f}')
    return rows


# ----- analysis ----------------------------------------------------------
def analyze(rows):
    """Return (onset, cross1pct, cross50pct, flags) for one system's rows."""
    onset = cross1 = cross50 = None
    for r in rows:
        if onset is None and r['bad'] > 0:
            onset = r['res']
        if cross1 is None and r['frac'] >= 0.01:
            cross1 = r['res']
        if cross50 is None and r['frac'] >= 0.50:
            cross50 = r['res']

    flags = []  # (res, kind, detail)
    for i, r in enumerate(rows):
        prev = rows[i - 1] if i > 0 else None
        # non-monotonic dip: DNC fell vs the immediately coarser resolution
        if prev is not None and r['frac'] + NOISE_TOL < prev['frac']:
            flags.append((r['res'], 'drop',
                          f'frac {100*r["frac"]:.2f}% < res {prev["res"]} '
                          f'{100*prev["frac"]:.2f}%'))
        # DNC island: DNC here but some finer resolution is clean
        if r['bad'] > 0 and any(rr['bad'] == 0 for rr in rows[i + 1:]):
            finer_clean = [rr['res'] for rr in rows[i + 1:] if rr['bad'] == 0]
            flags.append((r['res'], 'island',
                          f'DNC>0 but finer res {finer_clean} are clean'))
    return onset, cross1, cross50, flags


def report(all_rows):
    flagged_res = {}  # (sys,res) -> rows' dump, for cells we should write
    any_flag = False
    print('\n========== boundary + monotonicity report ==========')
    for sys, rows in all_rows.items():
        onset, c1, c50, flags = analyze(rows)
        print(f'\n{SYS_LABEL[sys]}: DNC onset at res {onset}; '
              f'>=1% at res {c1}; >=50% at res {c50}; '
              f'finest res {rows[-1]["res"]} = {100*rows[-1]["frac"]:.1f}% DNC')
        # always dump the onset resolution's cells for boundary characterization
        if onset is not None:
            flagged_res[(sys, onset)] = next(r for r in rows if r['res'] == onset)
        if flags:
            any_flag = True
            print(f'  !! {len(flags)} UNEXPECTED / non-monotonic flag(s):')
            for res, kind, detail in flags:
                print(f'     res {res} [{kind}] {detail}')
                flagged_res[(sys, res)] = next(r for r in rows if r['res'] == res)
        else:
            print('  monotonic: no unexpected DNC flags '
                  '(DNC only grows toward the finest resolutions).')

    # write flagged cells (onset + any flagged res) for reproduction
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with CELLS_FILE.open('w') as f:
        f.write('# system\tres\tcell_id\tgap\tvec3_vertices\n')
        for (sys, res), r in sorted(flagged_res.items()):
            for cid_str, gap, verts in r['dump']:
                vtxt = ';'.join(','.join(repr(x) for x in row) for row in verts) \
                    if verts is not None else 'NA'
                f.write(f'{sys}\t{res}\t{cid_str}\t{gap}\t{vtxt}\n')
    print(f'\nwrote flagged/onset DNC cells to {CELLS_FILE}')
    if not any_flag:
        print('RESULT: no unexpected DNCs — DNC is monotonic in resolution, '
              'consistent with the documented f64 floor (H3/ISEA7H stay clean; '
              'S2/A5 DNC only at the finest sub-metre levels).')
    else:
        print('RESULT: unexpected DNC pattern found — inspect the dump above.')


# ----- plot --------------------------------------------------------------
def plot(all_rows):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    for sys, rows in all_rows.items():
        res = [r['res'] for r in rows]
        frac = [100 * r['frac'] for r in rows]
        gap = [r['conv_worst'] if r['conv_worst'] > 0 else np.nan for r in rows]
        c = SYS_COLOR[sys]
        ax1.plot(res, frac, '-o', color=c, ms=4, label=SYS_LABEL[sys])
        ax2.plot(res, gap, '-o', color=c, ms=4, label=SYS_LABEL[sys])
    ax1.set_ylabel('DNC fraction (%)')
    ax1.set_title('DGGS did-not-converge fraction vs resolution '
                  '(skar default gap_tol = 1e-6)', fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax2.axhline(1e-6, color='0.4', ls='--', lw=1, label='gap_tol = 1e-6')
    ax2.set_yscale('log')
    ax2.set_ylabel('worst converged gap')
    ax2.set_xlabel('resolution / level')
    ax2.set_title('worst certified duality gap among converged cells '
                  '(f64 floor rises and crosses 1e-6 → DNC onset)', fontsize=11)
    ax2.grid(True, alpha=0.3, which='both')
    ax2.legend()
    fig.tight_layout()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG, dpi=200)
    plt.close(fig)
    print(f'wrote {PNG}')


def main():
    all_rows = {}
    for sys in ('h3', 's2', 'a5', 'isea7h'):
        t0 = time.perf_counter()
        all_rows[sys] = sweep_system(sys)
        total = sum(r['tested'] for r in all_rows[sys])
        print(f'[{sys}] {total:,} cells in {time.perf_counter() - t0:.1f}s')
    report(all_rows)
    plot(all_rows)


if __name__ == '__main__':
    main()

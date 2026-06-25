"""Cross-resolution DNC sweep — boundary + non-monotonicity, all systems.

Solves the cells of every resolution of H3, S2, A5, and the DGGAL ISEA7H/IVEA7H
hex grids with skar's DEFAULT settings (`skar.solve(v, geo='vec3')`,
gap_tol=1e-6) and characterizes where `did_not_converge` (DNC) appears:
  1. the DNC boundary (onset resolution and how the fraction rises);
  2. NON-MONOTONIC behaviour — a coarser resolution with more DNC than a finer
     one, or a "DNC island" (the H3-r7-r10-style surprise that v0.2.0 fixed);
  3. any unexpected DNCs, dumped reproducibly to out/dnc_sweep_cells.txt.

H3/ISEA7H/IVEA7H stay clean across all resolutions; S2/A5 DNC at their finest,
sub-metre resolutions (a genuine f64 duality-gap floor, ~22% of S2 L30 / ~47% of
A5 r30 at 1e-6 — see h3_gap_floor_report.md).

Reads the pre-generated cell sets (scripts/dggs/cells/, `just gen-cells`
first), so it imports no DGGS library and runs natively — every resolution at
up to N cells. Note it does not special-case the 12 H3 pentagons (a random set
rarely hits them); H3 is clean regardless.

Writes a 2-panel PNG (DNC fraction + worst converged gap vs resolution) to
out/dnc_sweep.png and prints per-system tables + a monotonicity report.

Run with:  just dnc-sweep
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

# _common.py (the cell-set reader) lives in the cells/ subfolder.
sys.path.insert(0, str(Path(__file__).resolve().parent / 'cells'))
import _common as cells  # noqa: E402

# ----- knobs -------------------------------------------------------------
NOISE_TOL = 1e-2          # monotonicity: ignore DNC-fraction dips below this
                          # (DNC lives in the finest, beyond-target resolutions,
                          # which are N_SMALL=25k -> ~0.3% sampling noise)
MAX_DUMP_PER_RES = 50     # cells written per flagged resolution

OUT_DIR = Path(__file__).resolve().parent / 'out'
PNG = OUT_DIR / 'dnc_sweep.png'
CELLS_FILE = OUT_DIR / 'dnc_sweep_cells.txt'

SYS_COLOR = {'h3': 'C0', 's2': 'C1', 'a5': 'C2', 'isea7h': 'C3', 'ivea7h': 'C4'}
SYS_LABEL = {'h3': 'H3', 's2': 'S2', 'a5': 'A5',
             'isea7h': 'ISEA7H', 'ivea7h': 'IVEA7H'}
SYSTEMS = list(SYS_COLOR)
# -------------------------------------------------------------------------


# ----- sweep -------------------------------------------------------------
def sweep_system(sys):
    """Return a list of per-resolution record dicts for `sys`'s cell sets."""
    rows = []
    print(f'\n=== {SYS_LABEL[sys]} ===')
    print(f'{"res":>3} {"tested":>10} {"dnc":>8} {"dnc%":>7} '
          f'{"conv_gap":>10} {"dnc_gap":>10} {"max_it":>6} {"secs":>6}')
    for res in cells.available_resolutions(sys):
        t0 = time.perf_counter()
        tested = dnc = infeas = raised = 0
        conv_worst = 0.0
        max_it = 0
        dnc_gaps = []
        dump = []  # (cid_str, gap, verts) for reproduction
        for cid, latlng in cells.load_cells(sys, res):
            tested += 1
            v = skar.to_vec3(latlng, geo='latlng_deg')
            try:
                r = skar.solve(v, geo='vec3')  # DEFAULT settings
            except Exception as e:
                raised += 1
                if len(dump) < MAX_DUMP_PER_RES:
                    dump.append((cid, f'raised:{type(e).__name__}', None))
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
                    dump.append((cid, r.gap, v))
        dt = time.perf_counter() - t0
        bad = dnc + infeas + raised
        frac = bad / tested if tested else 0.0
        dnc_gap_med = float(np.median(dnc_gaps)) if dnc_gaps else float('nan')
        rows.append(dict(sys=sys, res=res, tested=tested, dnc=dnc,
                         infeas=infeas, raised=raised, bad=bad, frac=frac,
                         conv_worst=conv_worst, dnc_gap_med=dnc_gap_med,
                         max_it=max_it, dump=dump))
        print(f'{res:>3} {tested:>10,} {bad:>8} {100*frac:>6.2f}% '
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
        # DNC island: a meaningful DNC fraction here but some finer resolution
        # is clean. Gated by NOISE_TOL so a stray cell or two at the f64 floor
        # (1 in 100k at sub-metre resolutions) isn't mistaken for a real band.
        if r['frac'] >= NOISE_TOL and any(rr['bad'] == 0 for rr in rows[i + 1:]):
            finer_clean = [rr['res'] for rr in rows[i + 1:] if rr['bad'] == 0]
            flags.append((r['res'], 'island',
                          f'frac {100*r["frac"]:.2f}% but finer res {finer_clean} clean'))
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
              'consistent with the documented f64 floor (H3/ISEA7H/IVEA7H stay '
              'clean; S2/A5 DNC only at the finest sub-metre levels).')
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
    for sys in SYSTEMS:
        t0 = time.perf_counter()
        all_rows[sys] = sweep_system(sys)
        total = sum(r['tested'] for r in all_rows[sys])
        print(f'[{sys}] {total:,} cells in {time.perf_counter() - t0:.1f}s')
    report(all_rows)
    plot(all_rows)


if __name__ == '__main__':
    main()

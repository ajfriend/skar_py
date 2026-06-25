"""DNC invariant check across every DGGS and resolution.

One pass/fail check (replacing dnc_sweep + dnc_stress): read the cached cell
sets, solve every cell with the DEFAULT solver settings, and verify the two
invariants that matter:

  1. clean where it's used — 0 DNC at each system's working (target) resolution
     and all coarser, i.e. the resolutions the survey and analyses actually use;
  2. monotone — DNC only grows toward the finest resolutions: no resolution has
     a meaningful DNC fraction while a finer one is clean (the H3-r7-r10-style
     "island" that skar_zig v0.2.0 removed), and the fraction never drops as
     resolution rises.

So `did_not_converge` appears only at the finest, sub-metre resolutions (S2/A5
onset ~r28) — the real f64 duality-gap floor — never at a working resolution and
never non-monotonically. Exits non-zero (printing the offenders) if either
invariant breaks.

(Unlike the old H3-only stress, this can't target the 12 pentagons specifically
— it tests whatever the cache sampled plus the fully-enumerated coarse
resolutions; H3 stays clean throughout regardless.)

Writes one diagnostic plot (out/dnc_check.png: DNC % vs resolution per system).
Reads the cell sets (`just gen-cells` first); native, no DGGS libraries.
Run with:  just dnc-check
"""

import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import skar

# _common.py (the cell-set reader) lives in the cells/ subfolder.
sys.path.insert(0, str(Path(__file__).resolve().parent / 'cells'))
import _common as cells  # noqa: E402

# ----- knobs -------------------------------------------------------------
NOISE_TOL = 1e-2          # DNC-fraction noise floor — the finest resolutions are
                          # N_SMALL=25k cells (~0.3% sampling noise), so a stray
                          # cell or two at the f64 floor isn't a real band
MAX_EXAMPLES = 5          # offending cell ids to print per failing resolution
SYS_COLOR = cells.SYS_COLOR
OUT = Path(__file__).resolve().parent / 'out' / 'dnc_check.png'
# -------------------------------------------------------------------------


def sweep_system(name):
    """[(res, tested, dnc, [example cids])] over the system's cached resolutions."""
    rows = []
    for res in cells.available_resolutions(name):
        tested = dnc = 0
        examples = []
        for cid, latlng in cells.load_cells(name, res):
            tested += 1
            r = skar.solve(skar.to_vec3(latlng, geo='latlng_deg'), geo='vec3')
            if not isinstance(r, skar.Converged):
                dnc += 1
                if len(examples) < MAX_EXAMPLES:
                    examples.append(cid)
        rows.append((res, tested, dnc, examples))
    return rows


def check_system(name, rows):
    """Return (failures, onset_res, finest_frac) for one system's sweep rows."""
    target = cells.TARGET_RES[name]
    frac = {res: dnc / tested for res, tested, dnc, _ in rows}
    reslist = [res for res, *_ in rows]
    failures = []

    for i, (res, tested, dnc, ex) in enumerate(rows):
        # invariant 1: clean at the working resolutions (res <= target)
        if res <= target and dnc:
            failures.append(f'r{res}: {dnc}/{tested} DNC at a working resolution '
                            f'(<= target r{target}); e.g. {ex}')
        # invariant 2a: a meaningful DNC fraction with a clean finer resolution
        if frac[res] >= NOISE_TOL and any(frac[r] == 0 for r in reslist[i + 1:]):
            clean = [r for r in reslist[i + 1:] if frac[r] == 0]
            failures.append(f'r{res}: {100*frac[res]:.1f}% DNC but finer res '
                            f'{clean} clean (non-monotone island)')
        # invariant 2b: fraction drops meaningfully going one step finer
        if i + 1 < len(rows):
            nxt = reslist[i + 1]
            if frac[nxt] + NOISE_TOL < frac[res]:
                failures.append(f'r{nxt}: {100*frac[nxt]:.1f}% DNC < r{res} '
                                f'{100*frac[res]:.1f}% (non-monotone drop)')

    onset = next((res for res, _, dnc, _ in rows if dnc), None)
    finest_res, finest_tested, finest_dnc, _ = rows[-1]
    return failures, onset, finest_dnc / finest_tested


def plot(per_system):
    """One diagnostic: DNC % vs resolution, a line per system. Flat-zero until
    the finest resolutions then a monotone rise — both invariants at a glance;
    an onset shift (the one regression the pass/fail check can't see) is visible
    here."""
    fig, ax = plt.subplots(figsize=(9, 5))
    for name, rows in per_system.items():
        res = [r for r, *_ in rows]
        pct = [100 * dnc / tested for _, tested, dnc, _ in rows]
        ax.plot(res, pct, '-o', ms=4, color=SYS_COLOR.get(name), label=name)
    ax.set_xlabel('resolution / level')
    ax.set_ylabel('did-not-converge (%)')
    ax.set_title('DGGS DNC fraction vs resolution (skar default gap_tol = 1e-6)')
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150)
    plt.close(fig)
    print(f'wrote {OUT}')


def main():
    print(f'{"system":8} {"onset":>6} {"finest %DNC":>12} {"result":>8}')
    all_failures = []
    per_system = {}
    for name in cells.TARGET_RES:
        t0 = time.perf_counter()
        rows = sweep_system(name)
        per_system[name] = rows
        failures, onset, finest_frac = check_system(name, rows)
        all_failures += [(name, f) for f in failures]
        print(f'{name:8} {("r%d" % onset if onset is not None else "none"):>6} '
              f'{100*finest_frac:>11.1f}% {"FAIL" if failures else "OK":>8} '
              f'  ({time.perf_counter() - t0:.0f}s)')

    plot(per_system)
    if all_failures:
        print('\nFAIL — DNC invariant(s) broken:')
        for name, f in all_failures:
            print(f'  {name:8} {f}')
        sys.exit(1)
    print('\nPASS: working resolutions clean; DNC only at the finest, '
          'sub-metre resolutions, and monotone.')


if __name__ == '__main__':
    main()

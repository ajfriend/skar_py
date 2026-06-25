"""Assert the solver converges at the working resolutions of every DGGS.

For each system, solves every cell of the *small* set at its target (working)
resolution and all coarser resolutions, and asserts none return
`did_not_converge` under the DEFAULT solver settings. This guards the invariant
that actually matters: the resolutions the survey and analyses use are
solver-safe.

(The finest S2/A5 resolutions *do* DNC — a real f64 duality-gap floor,
characterized by dnc_sweep.py — so this gate deliberately stops at each system's
working resolution, where everything converges with margin.)

Reads the pre-generated small cell sets (`just gen-cells` first); native
(skar + the cells group), no DGGS libraries, no Rosetta.

Run with:  just dggs-stress
Exits non-zero if any working-resolution cell fails to converge.
"""

import sys
import time
from pathlib import Path

import skar

# _common.py (the cell-set reader) lives in the cells/ subfolder.
sys.path.insert(0, str(Path(__file__).resolve().parent / 'cells'))
import _common as cells  # noqa: E402

# ----- knobs -------------------------------------------------------------
N_SMALL = 25_000          # small-set N the generators wrote
SEED = 0xC0FFEE
# Working (target) resolution per system — the finest in actual use. Cells at
# this resolution and all coarser must converge. (Mirrors survey.py / the
# generators' TARGET_RES.)
TARGET = {'h3': 9, 's2': 15, 'a5': 14, 'isea7h': 10, 'ivea7h': 10}
# -------------------------------------------------------------------------


def stress_system(name, target):
    """Solve every small-set cell at resolutions 0..target; (n_total, failures)."""
    total = 0
    bad = []  # (res, cid, status)
    for res in range(target + 1):
        for cid, latlng in cells.load_cells(name, res, N_SMALL, SEED):
            total += 1
            r = skar.solve(skar.to_vec3(latlng, geo='latlng_deg'), geo='vec3')
            if not isinstance(r, skar.Converged):
                bad.append((res, cid, getattr(r, 'status', 'error')))
    return total, bad


def main():
    print(f'{"sys":8} {"res<=":>5} {"cells":>11} {"bad":>6} {"secs":>6}')
    failures = 0
    for name, target in TARGET.items():
        t0 = time.perf_counter()
        total, bad = stress_system(name, target)
        failures += len(bad)
        print(f'{name:8} {target:>5} {total:>11,} {len(bad):>6} '
              f'{time.perf_counter() - t0:>6.1f}')
        for res, cid, status in bad[:10]:
            print(f'    !! r{res} {cid} -> {status}')

    if failures:
        print(f'\nFAIL: {failures} working-resolution cell(s) did not converge')
        sys.exit(1)
    print('\nPASS: every system converges at its working resolutions and coarser')


if __name__ == '__main__':
    main()

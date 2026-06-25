"""H3 DNC stress test across all resolutions, default solver settings.

Solves millions of H3 cells spanning every resolution (0..15) with
skar's DEFAULT settings — `skar.solve(verts, geo='vec3')`, no gap_tol or
max_outer overrides — and checks that none return `did_not_converge`
(or `infeasible`, or raise).

Per resolution:
  - always test all 12 pentagons (the most distorted cells, the likeliest
    to stress the solver);
  - if the resolution has few enough cells (<= ENUMERATE_MAX) enumerate
    EVERY cell (exhaustive at coarse resolutions, where random sampling
    would just hit the same handful repeatedly);
  - otherwise sample N_PER_RES uniform-random cells.

Any non-converged or erroring cell is recorded (id, status, gap) and its
id appended to OUT_FILE so the failure is reproducible. The script also
reports, per resolution, the worst achieved duality gap and max outer
iterations among the converged cells — i.e. how much headroom there is
below the 1e-6 default.

Total cells tested ~= sum(enumerated resolutions) + N_PER_RES * (number of
sampled resolutions). With the defaults below that's a few million; raise
N_PER_RES for a larger run.

Run with:  just dggs-stress
       (or: uv run --group dggs scripts/dggs/dnc_stress.py)
No CLI args (project convention) — edit the constants below in place.
"""

import time
from pathlib import Path

import numpy as np

import h3

import skar

# ----- knobs -------------------------------------------------------------
SEED = 0xC0FFEE
N_PER_RES = 500_000      # random cells per sampled (high) resolution
ENUMERATE_MAX = 300_000  # resolutions with <= this many cells: test all
CHUNK = 50_000           # sampling batch size (keeps memory flat)
RESOLUTIONS = range(0, 16)  # h3 supports 0..15
OUT_FILE = Path(__file__).resolve().parent / 'out' / 'dnc_failures.txt'
# -------------------------------------------------------------------------


def random_cells(res, n, rng):
    """Yield n uniform-on-sphere random cells at `res`, in batches."""
    done = 0
    while done < n:
        k = min(CHUNK, n - done)
        lng = 360.0 * rng.random(k) - 180.0
        lat = np.degrees(np.arcsin(2.0 * rng.random(k) - 1.0))  # equal-area
        for la, lo in zip(lat, lng):
            yield h3.latlng_to_cell(float(la), float(lo), res)
        done += k


def all_cells(res):
    """Yield every cell at `res` (via res-0 cells -> children)."""
    for c0 in h3.get_res0_cells():
        yield from h3.cell_to_children(c0, res)


def cells_for_res(res, rng):
    """The cell stream to test at `res`: pentagons + exhaustive-or-sampled."""
    yield from h3.get_pentagons(res)  # always, even when sampling
    if h3.get_num_cells(res) <= ENUMERATE_MAX:
        yield from all_cells(res)
    else:
        yield from random_cells(res, N_PER_RES, rng)


def solve_cell(cid):
    """(status, gap, outer_iters) for a cell under DEFAULT solver settings."""
    v = skar.to_vec3(h3.cell_to_boundary(cid), geo='latlng')
    r = skar.solve(v, geo='vec3')  # defaults: gap_tol=1e-6, max_outer=100
    return r


def main():
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    grand_total = grand_bad = 0
    failures = []  # (res, cid, status, gap)

    print(f'{"res":>3} {"mode":>10} {"tested":>10} {"bad":>6} '
          f'{"worst_gap":>10} {"max_it":>6} {"secs":>7}')
    for res in RESOLUTIONS:
        mode = 'all' if h3.get_num_cells(res) <= ENUMERATE_MAX else 'sample'
        t0 = time.perf_counter()
        tested = bad = 0
        worst_gap = 0.0
        max_it = 0
        for cid in cells_for_res(res, rng):
            tested += 1
            try:
                r = solve_cell(cid)
            except Exception as e:  # degenerate input, etc. — record, continue
                bad += 1
                failures.append((res, cid, f'raised:{type(e).__name__}', float('nan')))
                continue
            if isinstance(r, skar.Converged):
                worst_gap = max(worst_gap, r.gap)
                max_it = max(max_it, r.outer_iters)
            else:
                bad += 1
                failures.append((res, cid, r.status, getattr(r, 'gap', float('nan'))))
        dt = time.perf_counter() - t0
        grand_total += tested
        grand_bad += bad
        print(f'{res:>3} {mode:>10} {tested:>10,} {bad:>6} '
              f'{worst_gap:>10.2e} {max_it:>6} {dt:>7.1f}')

    print(f'\nTOTAL cells tested: {grand_total:,}')
    if grand_bad == 0:
        print('PASS — every cell converged under the default solver settings.')
        if OUT_FILE.exists():
            OUT_FILE.unlink()  # clear any stale failure log from a prior run
    else:
        print(f'FAIL — {grand_bad:,} cell(s) did not converge:')
        for res, cid, status, gap in failures[:20]:
            print(f'  r{res} {cid} {status} gap={gap:.3e}')
        if len(failures) > 20:
            print(f'  ... and {len(failures) - 20:,} more')
        with OUT_FILE.open('w') as f:
            for res, cid, status, gap in failures:
                f.write(f'{res}\t{cid}\t{status}\t{gap}\n')
        print(f'wrote all {grand_bad:,} failing cell ids to {OUT_FILE}')


if __name__ == '__main__':
    main()

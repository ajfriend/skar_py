"""Validate the corners-only enclosing-cone metric for the DGGAL grids.

The survey/sweep feed skar only a cell's *corner* vertices. For geodesic edges
that's exact: gnomonic maps geodesics to straight chords and the min-enclosing
ellipse of the corners already contains them (convexity). The equal-area DGGAL
grids have slightly *non-geodesic* edges that could bow outward at coarse
levels, so this checks it empirically for every grid in
`dggal_common.DGGAL_SYSTEMS`: across resolutions — including the coarsest levels
and the 12 pentagons (level 0) — it compares the aspect ratio from corners
against the ratio from edge-refined vertices (`getZoneRefinedWGS84Vertices`,
which densifies each edge). If the max delta is within solver tolerance,
corners-only is confirmed and we keep it.

Prints a report; writes nothing. Run under the dggs env:
    UV_PROJECT_ENVIRONMENT=.venv-dggs uv run --no-sync \
        scripts/dggs/validate_corners.py
"""

import numpy as np

import skar

import dggal_common

# ----- knobs -------------------------------------------------------------
SEED = 0xC0FFEE
REFINE = 20                 # edge-refinement points per edge for the reference
LEVELS = [0, 1, 2, 3, 5, 8, 11]   # incl coarsest + the 12 pentagons (level 0)
K = 300                     # cells tested per level (enumerate if fewer exist)
# -------------------------------------------------------------------------


def refined_vec3(ad, zone, refine):
    pts = ad.dggrs.getZoneRefinedWGS84Vertices(zone, refine)
    return skar.to_vec3(dggal_common.latlng_ring(pts), geo='latlng_deg')


def ar(verts):
    r = skar.solve(verts, geo='vec3')
    return r.aspect_ratio if isinstance(r, skar.Converged) else None


def cells_for_level(ad, level, rng):
    if ad.count(level) <= K:
        return list(ad.enumerate(level))
    seen, out = set(), []
    for zone in ad.sample(level, K * 4, rng):
        if zone not in seen:
            seen.add(zone)
            out.append(zone)
            if len(out) >= K:
                break
    return out


def check(name, ad):
    """Report corners-vs-refined max |dAR| per level; return the overall max."""
    rng = np.random.default_rng(SEED)
    print(f'\n{name} corners-vs-refined (edgeRefinement={REFINE})')
    print(f'{"lvl":>3} {"cells":>6} {"pents":>6} {"max|dAR|":>10} '
          f'{"max_rel":>10} {"corners_AR_range":>22}')
    overall = 0.0
    for level in LEVELS:
        zones = cells_for_level(ad, level, rng)
        npent = 0
        max_abs = max_rel = 0.0
        ars = []
        for z in zones:
            corners = ad.verts(z)
            if corners.shape[0] == 5:
                npent += 1
            a_c = ar(corners)
            a_r = ar(refined_vec3(ad, z, REFINE))
            if a_c is None or a_r is None:
                continue
            ars.append(a_c)
            d = abs(a_c - a_r)
            max_abs = max(max_abs, d)
            max_rel = max(max_rel, d / a_r)
        overall = max(overall, max_abs)
        rng_txt = f'[{min(ars):.4f}, {max(ars):.4f}]' if ars else 'n/a'
        print(f'{level:>3} {len(zones):>6} {npent:>6} {max_abs:>10.2e} '
              f'{max_rel:>10.2e} {rng_txt:>22}')
    return overall


def main():
    worst = 0.0
    for _s in dggal_common.DGGAL_SYSTEMS.values():
        worst = max(worst, check(_s['cls'], dggal_common.Adapter(_s['cls'])))
    print(f'\noverall max |dAR| across all grids = {worst:.3e}')
    print('corners-only CONFIRMED (within solver tolerance)' if worst < 1e-3
          else 'corners-only delta NOT negligible — investigate')


if __name__ == '__main__':
    main()

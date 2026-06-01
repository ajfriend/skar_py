"""Country aspect-ratio example.

Loads the Natural Earth admin-0 countries with geopandas, computes each
country's tightest enclosing-cone aspect ratio with `skar`, prints the full
ranking, and writes one PNG per country (boundary + enclosing ellipse) into
out/<slug>.png.

Like the states example, it's one pass: geopandas hands us shapely geometries,
`skar.solve(geom)` consumes them via `__geo_interface__`, and `skar.plot_cone`
draws the outline + enclosing ellipse — no vertex extraction, no intermediate
files. Countries can be harder than US states — France's territory reaches
across to French Guiana, and very elongated countries like Chile (its mainland
spans ~38° of latitude) are slow to certify — so a few need more than skar's
default 100 outer iterations to certify the gap; we raise `max_outer`. Anything
that still doesn't converge, is proven `Infeasible` (no hemisphere holds it),
or is too degenerate to solve (ValueError) is reported and skipped — not a
failure.

Run with:  just countries  (or: uv run --group geo scripts/countries/countries.py)
No CLI args (project convention) — edit the constants below in place.
"""

import re
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import geopandas as gpd

import skar

# ----- knobs -------------------------------------------------------------
URL = ('https://raw.githubusercontent.com/nvkelso/natural-earth-vector/'
       'master/geojson/ne_110m_admin_0_countries.geojson')
NAME_FIELD = 'ADMIN'        # canonical, always-present country name
MAX_OUTER = 1000            # a few elongated countries need >100 to certify
OUT_DIR = Path(__file__).resolve().parent / 'out'
DPI = 200
# -------------------------------------------------------------------------


def slug(name):
    return re.sub(r'[^a-z0-9]+', '_', name.lower()).strip('_')


def save_plot(name, geom, r):
    ax = skar.plot_cone(r, geom, title=f'{name} — AR {r.aspect_ratio:.3f}')
    ax.figure.savefig(OUT_DIR / f'{slug(name)}.png', dpi=DPI, bbox_inches='tight')
    plt.close(ax.figure)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    gdf = gpd.read_file(URL)

    converged = []
    for name, geom in zip(gdf[NAME_FIELD], gdf.geometry):
        try:
            r = skar.solve(geom, max_outer=MAX_OUTER)
        except ValueError as e:                        # degenerate / too few pts
            print(f'  skip {name}: {e}')
            continue
        if isinstance(r, skar.Converged):
            converged.append((str(name), geom, r))
        else:
            print(f'  skip {name}: {r.status}')
    converged.sort(key=lambda t: t[2].aspect_ratio, reverse=True)

    print(f'\n{"rank":>4}  {"AR":>7}  country   (most to least elongated)')
    for i, (name, _, r) in enumerate(converged, 1):
        print(f'{i:>4}  {r.aspect_ratio:>7.3f}  {name}')

    for name, geom, r in converged:
        save_plot(name, geom, r)
    print(f'\nwrote {len(converged)} plots to {OUT_DIR}')


if __name__ == '__main__':
    main()

"""US-state aspect-ratio example.

Loads the US-states GeoJSON with geopandas, computes each state's tightest
enclosing-cone aspect ratio with `skar`, prints the full ranking, and writes
one PNG per state (boundary + enclosing ellipse) into out/<slug>.png.

The whole thing is one pass: geopandas hands us shapely geometries, and
`skar.solve(geom)` consumes them directly via `__geo_interface__` (the
MultiPolygon's exterior rings become the point set). `skar.plot_cone` then
draws the outline + enclosing ellipse — no vertex extraction, no intermediate
files. The old gen -> JSON -> Zig -> JSON -> plot pipeline collapses to
load -> solve -> plot.

Run with:  just states     (or: uv run --group geo scripts/states/states.py)
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
URL = ('https://raw.githubusercontent.com/PublicaMundi/MappingAPI/'
       'master/data/geojson/us-states.json')
EXCLUDE = {'District of Columbia', 'Puerto Rico'}  # land on exactly 50
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
    gdf = gdf[~gdf['name'].isin(EXCLUDE)]

    converged = []
    for name, geom in zip(gdf['name'], gdf.geometry):
        r = skar.solve(geom, max_outer=1000)
        if isinstance(r, skar.Converged):
            converged.append((name, geom, r))
        else:
            print(f'  skip {name}: {r.status}')
    converged.sort(key=lambda t: t[2].aspect_ratio, reverse=True)

    print(f'\n{"rank":>4}  {"AR":>7}  state   (most to least elongated)')
    for i, (name, _, r) in enumerate(converged, 1):
        print(f'{i:>4}  {r.aspect_ratio:>7.3f}  {name}')

    for name, geom, r in converged:
        save_plot(name, geom, r)
    print(f'\nwrote {len(converged)} plots to {OUT_DIR}')


if __name__ == '__main__':
    main()

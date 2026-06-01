"""Country aspect-ratio example.

Loads the Natural Earth admin-0 countries with geopandas, computes each
country's tightest enclosing-cone aspect ratio with `skar`, prints the full
ranking, and writes one PNG per country (boundary + enclosing ellipse) into
out/<slug>.png.

Like the states example, it's one pass: geopandas hands us shapely geometries
and `skar.solve(geom)` consumes them via `__geo_interface__` — no vertex
extraction, no intermediate files. Countries are larger and more transoceanic
than US states (France with French Guiana, Chile with Easter Island), so a few
need more than skar's default 100 outer iterations to certify the gap; we raise
`max_outer`. Anything that still doesn't converge, is proven `Infeasible` (no
hemisphere holds it), or is too degenerate to solve (ValueError) is reported
and skipped — not a failure.

Run with:  just countries  (or: uv run --group geo scripts/countries/countries.py)
No CLI args (project convention) — edit the constants below in place.
"""

import re
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import geopandas as gpd

import skar

# ----- knobs -------------------------------------------------------------
URL = ('https://raw.githubusercontent.com/nvkelso/natural-earth-vector/'
       'master/geojson/ne_110m_admin_0_countries.geojson')
NAME_FIELD = 'ADMIN'        # canonical, always-present country name
MAX_OUTER = 1000            # a few elongated countries need >100 to certify
EARTH_R_M = 6_371_008.8
OUT_DIR = Path(__file__).resolve().parent / 'out'
DPI = 200
# -------------------------------------------------------------------------


def slug(name):
    return re.sub(r'[^a-z0-9]+', '_', name.lower()).strip('_')


def rings_lonlat(geom):
    """Per-ring (lon, lat) arrays for a shapely (Multi)Polygon — exterior and
    holes, each ring kept separate so the plot doesn't draw spurious segments
    across disjoint pieces (Russia, Indonesia, Canada, ...)."""
    polys = geom.geoms if geom.geom_type == 'MultiPolygon' else [geom]
    return [np.asarray(ring.coords)[:, :2]
            for poly in polys for ring in [poly.exterior, *poly.interiors]]


def save_plot(name, geom, r):
    """Gnomonic-project the outline at the cone axis, overlay the enclosing
    ellipse, and write out/<slug>.png, in skar's eigenbasis Q (major axis
    horizontal, semi-axes sqrt(2/3)/sigma[1:])."""
    b, U = r.Q[:, 0], r.Q[:, 1:]
    semi = np.sqrt(1.0 - r.sigma[0] ** 2) / r.sigma[1:] * EARTH_R_M
    # Orient north-up: a 180° turn (keeps the major axis horizontal) if world
    # north would otherwise project downward. Q is right-handed, so the
    # outline's chirality is already correct — only this flip is ambiguous.
    north = np.array([0.0, 0.0, 1.0]) - b[2] * b
    flip = -1.0 if (north @ U)[1] < 0.0 else 1.0

    fig, ax = plt.subplots(figsize=(7, 7))
    first = True
    for ring in rings_lonlat(geom):
        v = skar.to_vec3(ring, geo='lonlat')           # GeoJSON order, no swap
        y = (v @ U) / (v @ b)[:, None] * flip * EARTH_R_M
        closed = np.vstack([y, y[:1]])
        ax.plot(closed[:, 0], closed[:, 1], '-', color='C0', lw=0.9,
                label='boundary' if first else None)
        first = False

    t = np.linspace(0.0, 2.0 * np.pi, 400)
    ax.plot(semi[0] * np.cos(t), semi[1] * np.sin(t), '-', color='0.25', lw=1.5,
            label='enclosing ellipse')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('major axis (m)')
    ax.set_ylabel('minor axis (m)')
    ax.text(0.03, 0.97, f'{name}\nAR {r.aspect_ratio:.4f}', transform=ax.transAxes,
            va='top', ha='left', fontsize=10,
            bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.85))
    ax.legend(loc='lower right', fontsize=8)
    fig.suptitle(f'{name}: tightest enclosing cone (AR {r.aspect_ratio:.3f})')
    fig.tight_layout()
    fig.savefig(OUT_DIR / f'{slug(name)}.png', dpi=DPI)
    plt.close(fig)


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

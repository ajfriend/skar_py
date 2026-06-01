# Country aspect ratios

Compute the tightest enclosing-cone aspect ratio for each country and plot the
most-elongated ones.

```sh
just countries
# or: uv run --group geo scripts/countries/countries.py
```

Prints the full ranking (most → least elongated) and writes one PNG per country
to `scripts/countries/out/<slug>.png` (gitignored) — each country
gnomonic-projected at its cone axis with the enclosing ellipse overlaid. (France
and Chile top the ranking — French Guiana and Easter Island stretch their
enclosing cones.)

## Design

Same single-pass shape as the states example: `geopandas` loads the Natural
Earth admin-0 countries and `skar.solve(geom)` consumes each shapely geometry
via `__geo_interface__`. Differences:

- `max_outer = 1000` — a few transoceanic countries (France, Chile) need more
  than skar's default 100 outer iterations to certify the gap at strict
  tolerance; the aspect ratio is stable well before then.
- Outcomes that are `Infeasible` (no hemisphere holds the country) or too
  degenerate to solve (`ValueError`) are reported and skipped, not fatal.

Config (`URL`, `NAME_FIELD`, `MAX_OUTER`) is at the top of `countries.py`.

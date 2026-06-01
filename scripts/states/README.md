# US-state aspect ratios

Compute the tightest enclosing-cone aspect ratio for each US state and plot
the most-elongated ones.

```sh
just states
# or: uv run --group geo scripts/states/states.py
```

Prints the full ranking (most → least elongated) and writes one PNG per state
to `scripts/states/out/<slug>.png` (gitignored) — each state gnomonic-projected
at its cone axis with the enclosing ellipse overlaid, major axis horizontal.

## Design

One pass, no intermediate files. `geopandas` loads the US-states GeoJSON, and
each shapely geometry is fed straight to `skar.solve(geom)` — the geometry's
`__geo_interface__` turns its (multi)polygon rings into the point set. The plot
re-projects each ring with `skar.to_vec3(ring, geo='lonlat')` (GeoJSON axis
order). This collapses the original Zig-repo pipeline
(`gen_states.py → states.json → states.zig → states_aspect.json → plot`) into
a single script now that `skar.solve` is callable from Python.

Config (`EXCLUDE`, `URL`) is in constants at the top of `states.py`.

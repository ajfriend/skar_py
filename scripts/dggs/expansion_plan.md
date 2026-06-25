# DGGS comparison expansion plan

Tracking issue: [#6](https://github.com/ajfriend/skar_py/issues/6).

Goal: extend the DGGS aspect-ratio comparison beyond **H3 / S2 / A5** to span
more cell shapes, projections, and apertures. Work through the target systems
**one at a time** — each is its own PR and leaves the scripts runnable.

This doc is the standing spec for that work. A fresh session should be able to
pick the next unchecked system and execute it end-to-end from here.

## Target systems

Coverage goal: **{ISEA, IVEA} projections × {3, 7} apertures** for hexagons,
plus a **quad** (rHEALPix) and a **triangle** (ISEA4T) cell shape, against the
existing H3 / S2 / A5 baselines.

| System   | Cell shape | Projection | Aperture | Binding                            | Done |
|----------|------------|------------|----------|------------------------------------|------|
| ISEA3H   | hexagon    | ISEA       | 3        | DGGAL `ISEA3H`                     | [ ]  |
| ISEA7H   | hexagon    | ISEA       | 7        | DGGAL `ISEA7H`                     | [x]  |
| IVEA3H   | hexagon    | IVEA       | 3        | DGGAL `IVEA3H`                     | [ ]  |
| IVEA7H   | hexagon    | IVEA       | 7        | DGGAL `IVEA7H`                     | [x]  |
| IGEO7    | hexagon    | ISEA       | 7        | DGGAL `ISEA7H_Z7`                  | [ ]  |
| rHEALPix | quad       | HEALPix    | 3        | DGGAL `RHEALPix`                   | [ ]  |
| ISEA4T   | triangle   | ISEA       | 4        | dggrid4py (portable DGGRID binary) | [ ]  |

**Note — IGEO7 vs ISEA7H geometry.** `ISEA7H_Z7` is ISEA7H with Z7 *addressing*;
the cells are geometrically identical. The aspect-ratio metric reads geometry,
not cell IDs, so IGEO7 and ISEA7H will produce the **same** distribution.
Decide before doing both: either keep one (label it IGEO7, the recognized name)
or keep both and document that they coincide. Don't silently ship two identical
curves.

## Binding strategy

Six of seven come from **one library**: DGGAL (`pip install dggal`). One DGGAL
helper module — `dggal_common.py`, written for ISEA7H — is reused for all six;
each new system is a one-liner `Adapter('<DGGRSClass>')`.

> **Platform caveat (discovered wiring ISEA7H).** DGGAL's prebuilt macOS
> **arm64** wheel (dggal/ecrt 0.0.6) is arch-broken: it bundles **x86_64**
> `libecrt.dylib`/`libdggal.dylib` inside an arm64 wheel, so it can't load on
> Apple Silicon. Correct arm64 wheels exist only for cp310 (0.0.4–0.0.5), and
> skar needs ≥3.11, so there is **no native arm64 + py≥3.11 combo**. Even the
> `ecdev` eC toolchain's arm64 wheel ships x86_64 binaries, so building from
> source means bootstrapping the whole Ecere SDK cross-arch (impractical). The
> fix is published only in the unreleased 0.0.7rc1. **Workaround in this repo:**
> the `dggs` just-targets run under an **x86_64 (Rosetta) Python 3.13** in a
> separate `.venv-dggs` (the x86_64 and Linux wheels are self-consistent); the
> native arm64 dev env is untouched. Revisit once a fixed arm64 wheel ships.

ISEA4T (triangles) is the only system DGGAL doesn't expose; it comes from
**dggrid4py** with `tool.get_portable_executable()` so no system DGGRID install
is required. Do ISEA4T last.

## Design decisions

- **Corners only — no edge densification.** The enclosing-cone metric is the
  minimum enclosing ellipse of the gnomonic image of the points; gnomonic maps
  geodesic edges to straight chords, and an ellipse enclosing the corners
  encloses those chords (convexity), so densified edge points are never active.
  This holds exactly for geodesic edges; the equal-area grids have
  non-geodesic edges that *could* bow outward at coarse resolutions. So:
  ship corners-only, but **validate once** (see first-system tasks).
- **Same as the existing adapters**: vertices are delivered as an `(M, 3)`
  array of **unit** vec3 (M varies — 6 for hexagons, 5 for the 12 pentagons,
  4 for quads, 3 for triangles), no repeated closing vertex.
- **No CLI args** (project convention): all knobs are constants edited in place.

## Where each system plugs in

DGGAL grids (every row except ISEA4T) are driven by a single registry,
`dggal_common.DGGAL_SYSTEMS`. **Adding one is a single row** + a calibrate run +
docs; calibrate.py / survey.py / dnc_sweep.py / validate_corners.py all loop the
registry, so there are no per-system functions to write.

### 1. `dggal_common.DGGAL_SYSTEMS` — the one code edit
Add a row:
```python
'<key>': dict(cls='<DGGRSClass>', color='C<n>', res=<level>, scan=range(...)),
```
`cls` = DGGAL class; `color` = next matplotlib slot; `res` = H3-r9-matched level
(fill in after step 2); `scan` = calibrate search range.

### 2. `just calibrate`
calibrate.py scans the new grid automatically (it loops the registry). Read off
the picked level and set it as the row's `res`.

### 3. `just dggs` / `just dnc-sweep`
survey.py (label/color/iterator) and dnc_sweep.py (the six-method adapter entry
+ `res_range` from `max_level()`, color, label, `N_PER_RES`) are derived from the
registry; `main()` iterates `SYSTEMS`. Nothing else to touch.

### 4. `validate_corners.py`
Loops the registry, so the new grid's corners-only check runs automatically.

### 5. Docs
`scripts/dggs/README.md` system list / coverage, the repo `readme.md` tree line,
and a terse `changelog.md` bullet.

> **pyproject.toml**: `dggal` is already in the `dggs` group — no change for
> DGGAL grids. **ISEA4T is the exception**: it isn't a DGGAL DGGRS, so it comes
> from `dggrid4py` (add that dep) with its own adapter, not a `DGGAL_SYSTEMS`
> row.

## Shared DGGAL helper — `dggal_common.py` (built)

Built during ISEA7H (#7) and extended for IVEA7H. `scripts/dggs/dggal_common.py`:

- initializes the DGGAL `Application` **once** (`pydggal_setup`) at import
  (with a guarded `dlopen` fallback for the broken arm64 wheel);
- `Adapter(cls)` wraps a DGGRS and exposes `count` / `enumerate` / `sample` /
  `verts` / `cid_str` / `max_level` / `area_km2` / `iter_sample`;
- `latlng_ring(points)` converts DGGAL WGS84 vertices to `(lat, lng)` (corners
  only, closing-repeat stripped) for `skar.to_vec3(..., geo='latlng_deg')`;
- `DGGAL_SYSTEMS` is the registry the four scripts loop over (see above).

Map onto the DGGAL `DGGRS` API:

| Need              | DGGAL call                                                       |
|-------------------|------------------------------------------------------------------|
| valid resolutions | `range(0, dggrs.getMaxDGGRSZoneLevel() + 1)` (cap for the sweep) |
| `count(res)`      | `dggrs.countZones(level)`                                        |
| `enumerate(res)`  | `dggrs.listZones(level, worldBBox)`                              |
| `verts(zone)`     | `dggrs.getZoneWGS84Vertices(zone)` → lat/lng → vec3              |
| `cid_str(zone)`   | `dggrs.getZoneTextID(zone)`                                      |

**Unknowns — all resolved wiring ISEA7H (#7):**
1. Whole-world bbox: `wholeWorld` global (`GeoExtent((-90,-180),(90,180))`).
2. Point→zone: `getZoneFromWGS84Centroid(level, GeoPoint(lat, lng))`.
3. Pentagons handled — `verts`/`iter_sample` read the returned vertex count, so
   5-vertex pentagons and 6-vertex hexagons both work.
4. dggal is **BSD-3-Clause** — compatible. (But its macOS arm64 wheel is
   arch-broken; see the platform caveat above.)

## Per-system recipe (DGGAL grids)

1. Branch from `main`.
2. Add one row to `dggal_common.DGGAL_SYSTEMS` (`cls`, `color`, placeholder
   `res`, `scan`).
3. `just calibrate` → set the row's `res` to the picked level.
4. `just dggs` → confirm the new curve appears and solves cleanly (0 DNC) at the
   matched resolution.
5. Standalone-sweep the new grid (and regress an existing one) — confirm 0 DNC /
   monotonic. Skip the ~45-min full combined sweep; its PNG is gitignored and
   H3/S2/A5 are unchanged.
6. `uv run … validate_corners.py` → corners-only confirmed (loops the registry).
7. Update the READMEs + a terse `changelog.md` bullet; tick the row above.
8. `just test`; open a PR (full detail in the PR, changelog points at it).

(ISEA4T is **not** DGGAL — it follows the older per-adapter recipe with
`dggrid4py`, so it adds the dep + its own adapter rather than a registry row.)

## First-system extra task — validate corners-only

Once the first DGGAL system runs, add a throwaway validation script
(`scripts/dggs/validate_corners.py`, gitignored output) that, for a sample of
cells across resolutions (include the **coarsest** levels and the **pentagons**),
compares the aspect ratio from corners vs. from
`getZoneRefinedWGS84Vertices(zone, n)`. If the max delta is within solver
tolerance, corners-only is confirmed empirically and we keep it everywhere.
Record the result (a sentence in the PR / this doc) and move on.

> **Result (ISEA7H + IVEA7H).** `validate_corners.py` compares corners vs.
> `getZoneRefinedWGS84Vertices(zone, 20)` across levels 0/1/2/3/5/8/11 (level 0
> = the 12 pentagons) for every registry grid: overall max ΔAR ≈ **2e-6**, i.e.
> at the 1e-6 gap floor, not systematic edge-bowing. Corners-only confirmed;
> kept. New DGGAL grids are validated automatically (the script loops the
> registry).

## Conventions (see also global CLAUDE.md + repo memory)

- `uv` only; run scripts with `uv run --group dggs ...` or the `just` targets.
- No CLI arg parsing — knobs are in-file constants.
- Single quotes for Python strings.
- Land changes via PRs; keep changelog/release notes terse.
- These systems are added on their own technical merit — do not attribute the
  selection to any external person or private communication.

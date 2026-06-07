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
| IVEA7H   | hexagon    | IVEA       | 7        | DGGAL `IVEA7H`                     | [ ]  |
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

Adding a DGGS touches four files. All three scripts hardcode the system list,
so each new system needs an entry in each.

### 1. `pyproject.toml` — `dggs` dependency group (~line 59)
Add the binding (`dggal`, and later `dggrid4py`) to the `dggs` group.

### 2. `scripts/dggs/calibrate.py` — area matching
Add `<sys>_area(res, n) -> median area in km^2` (use `sparea`, skar-free),
register it in `AREA_FN`, add a `SCAN` range. Re-run `just calibrate`, then
bake the printed pick into `survey.py`.

### 3. `scripts/dggs/survey.py` — the survey at a matched resolution
Add `iter_<sys>(n, seed)` yielding `(id_str, (M, 3) unit-vertex array)`,
register in `ITERATORS`, append to `SYSTEMS`, `SYS_LABEL`, `SYS_COLOR`, and add
the baked-in `<SYS>_RES` constant from calibrate.

### 4. `scripts/dggs/dnc_sweep.py` — full cross-resolution DNC sweep
Add the six-function adapter and register it in `SYSTEMS`:
```python
'<sys>': dict(
    count=...,      # count(res) -> int
    enumerate=...,  # enumerate(res) -> iterator of cell ids
    sample=...,     # sample(res, n, rng) -> iterator of cell ids
    verts=...,      # verts(cid) -> (M, 3) unit vec3, corners only
    cid_str=...,    # cid_str(cid) -> str
    res_range=...,  # range of valid resolutions
)
```
Also add to `SYS_COLOR`, `SYS_LABEL`, and the system tuple in `main()`.

### 5. `scripts/dggs/README.md` (and repo `readme.md` if it lists systems)
Update the system list / coverage summary.

## Shared DGGAL helper (build during the first DGGAL system)

Create `scripts/dggs/dggal_common.py` so the six DGGAL systems don't each
re-derive the binding glue. It should:

- initialize the DGGAL `Application` **once** (`pydggal_setup`) at import;
- given a DGGRS instance, expose: `count(level)`, `enumerate(level)`,
  `sample(level, n, rng)`, `verts(zone) -> (M,3) unit vec3`,
  `cid_str(zone)`, and `area(level, n)` for calibrate;
- convert `getZoneWGS84Vertices(zone)` (WGS84 lat/lon, corners only) to unit
  vec3 — reuse `skar.to_vec3(..., geo='latlng_deg')` with `(lat, lon)` tuples,
  matching how the A5 adapter remaps `(lon, lat)` rings today.

Map onto the DGGAL `DGGRS` API:

| Need              | DGGAL call                                                       |
|-------------------|------------------------------------------------------------------|
| valid resolutions | `range(0, dggrs.getMaxDGGRSZoneLevel() + 1)` (cap for the sweep) |
| `count(res)`      | `dggrs.countZones(level)`                                        |
| `enumerate(res)`  | `dggrs.listZones(level, worldBBox)`                              |
| `verts(zone)`     | `dggrs.getZoneWGS84Vertices(zone)` → lat/lon → vec3              |
| `cid_str(zone)`   | `dggrs.getZoneTextID(zone)`                                      |

**Unknowns to resolve while wiring the first system** (ISEA7H or ISEA3H):
1. Exact construction of the whole-world bbox for `listZones`.
2. The point→zone lookup for sampling (uniform lon/lat → containing zone);
   confirm the method name in the installed binding. Fallback: `listZones` at
   the level, then `rng`-sample the returned list (fine for coarse levels).
3. Pentagon handling — the 12 pentagons return 5 vertices; make sure `verts`
   and sampling don't assume 6.
4. Confirm DGGAL's license is compatible with skar before adding the dep.

## Per-system recipe (repeat for each row)

1. Branch from `main`.
2. Add the binding to the `dggs` group in `pyproject.toml` (first DGGAL
   system only; ISEA4T adds `dggrid4py`).
3. Wire the system into `calibrate.py`, run `just calibrate`, bake the picked
   resolution into `survey.py`.
4. Wire it into `survey.py` and `dnc_sweep.py`.
5. Run `just dggs` and `just dnc-sweep`; confirm the new system appears, solves
   cleanly at the matched resolution, and its DNC behaviour is monotonic (same
   bar the documented sub-metre f64 floor).
6. Update the READMEs.
7. Open a PR. Per project convention, the PR carries the full detail; keep the
   changelog terse and pointing at the PR/commit.

## First-system extra task — validate corners-only

Once the first DGGAL system runs, add a throwaway validation script
(`scripts/dggs/validate_corners.py`, gitignored output) that, for a sample of
cells across resolutions (include the **coarsest** levels and the **pentagons**),
compares the aspect ratio from corners vs. from
`getZoneRefinedWGS84Vertices(zone, n)`. If the max delta is within solver
tolerance, corners-only is confirmed empirically and we keep it everywhere.
Record the result (a sentence in the PR / this doc) and move on.

> **Result (ISEA7H).** `validate_corners.py` compared corners vs.
> `getZoneRefinedWGS84Vertices(zone, 20)` across levels 0/1/2/3/5/8/11
> (level 0 = the 12 pentagons): overall max ΔAR ≈ **1.6e-6**, i.e. at the
> 1e-6 gap floor, not systematic edge-bowing. Corners-only confirmed; kept.

## Conventions (see also global CLAUDE.md + repo memory)

- `uv` only; run scripts with `uv run --group dggs ...` or the `just` targets.
- No CLI arg parsing — knobs are in-file constants.
- Single quotes for Python strings.
- Land changes via PRs; keep changelog/release notes terse.
- These systems are added on their own technical merit — do not attribute the
  selection to any external person or private communication.

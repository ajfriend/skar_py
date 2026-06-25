# DGGS cell cache + analyses

The Parquet-based DGGS pipeline: **generate** random cells for each system into
Parquet once, then run every analysis **off those files** — natively, with no
DGGS library and no x86_64/Rosetta env. (The live, can't-be-cached analyses —
edge refinement, grid/neighbor explorations — live in `../dggs_old/`.)

## Layout

```
cells/            generation — one PEP 723 / uv-run script per DGGS + the cache core
  _common.py        schema, generate(), load_cells(), and the pipeline config
                    (SEED, N_BIG/N_SMALL, TARGET_RES)
  gen_{h3,s2,a5,dggal}.py   (gen_dggal inlines the bit of the DGGAL engine it
                            needs; the fuller glue is in ../dggs_old/)
  out/              the Parquet cell sets (gitignored): {dggs}_r{res}.parquet
survey.py         per-system aspect-ratio survey (-> out/histograms.png, extremes.png)
calibrate.py      match S2/A5/DGGAL resolutions to an H3 r9 cell by area
dnc_check.py      DNC invariants: working resolutions clean + DNC only at the
                  finest, monotone (one pass/fail check, every system)
explorations/     ar_histograms.py, ar_vs_pca.py (cache-reading)
```

## Use

```sh
just gen-cells     # generate all the Parquet (run once; ~2 min, dggal under Rosetta)
just dggs          # survey -> out/histograms.png + extremes.png
just calibrate     # area-match resolutions (pick -> bake into _common.TARGET_RES)
just dnc-check     # assert the DNC invariants (clean working res, monotone)
```

Each analysis is a plain `load_cells(dggs, res)` over the cache + `skar.solve`;
none import h3/s2sphere/a5_fast/dggal. See `cells/README.md` for the schema,
filename convention, and how a resolution's cell budget is chosen.

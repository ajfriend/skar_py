# Cache-reading explorations

Investigative AR plots that read the Parquet cell sets (no DGGS library, native):

- **`ar_histograms.py`** — ISEA7H vs IVEA7H aspect-ratio distributions at the
  H3-r9-matched level (ISEA7H is right-skewed with a long tail; IVEA7H is a
  flat-topped, sharply-bounded band). → `out/ar_histograms.png`
- **`ar_vs_pca.py`** — per-cell skar enclosing-cone AR vs PCA second-moment AR;
  shows the bulk lies on the `skar == PCA` diagonal (the metric is the real
  smooth distortion field) while the low-skar tail is grid accidents, not round
  cells. → `out/ar_vs_pca.png`

```sh
uv run --group cells scripts/dggs_cache/explorations/ar_histograms.py
```

`_common.py` provides the skar primitives (`aspect_ratio`, `gnomonic_xy`,
`tangent_basis_vec`, ...) plus a lazy `cells` handle on the cache reader
(`../cells/_common.py`). The grid/neighbor explorations that need the live DGGAL
engine are in `../../dggs_old/explorations/`.

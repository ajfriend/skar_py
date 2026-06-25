# /// script
# requires-python = "==3.13.*"
# dependencies = ["dggal", "numpy>=1.24", "pyarrow>=15"]
# ///
"""Generate random DGGAL cell sets (ISEA7H, IVEA7H, ...) -> Parquet.

One file per system in dggal_common.DGGAL_SYSTEMS per level 0..target (coarse
levels with <= N cells are enumerated in full, finer ones sampled; schema +
write logic in _common.py). Reuses the parent dggal_common.py for the DGGAL
setup + registry, so adding a system is still one row there.

dggal ships an arch-broken macOS arm64 wheel, so on Apple Silicon run under an
x86_64 (Rosetta) Python 3.13 — the wheel is self-consistent there:

    uv run --python cpython-3.13-macos-x86_64 scripts/dggs/cells/gen_dggal.py

On Linux (correct wheels) a plain `uv run` works. No CLI args (project
convention) — edit the constants below in place.
"""

import sys
from pathlib import Path

# dggal_common.py lives one level up (scripts/dggs/); reuse its DGGAL glue.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import dggal_common as dc  # noqa: E402

import _common  # noqa: E402

# ----- knobs -------------------------------------------------------------
N = 100_000
SEED = 0xC0FFEE
# -------------------------------------------------------------------------


if __name__ == '__main__':
    for name, row in dc.DGGAL_SYSTEMS.items():
        ad = dc.Adapter(row['cls'])
        _common.generate_levels(
            name, row['res'], N, SEED,
            zone_at=lambda res, lon, lat, _ad=ad: _ad.zone_at(res, lon, lat),
            cid_str=ad.cid_str,
            ring_of=ad.ring_latlng,
            count_at=lambda res, _ad=ad: _ad.count(res),
            enumerate_at=lambda res, _ad=ad: _ad.enumerate(res),
        )

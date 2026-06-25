# /// script
# requires-python = "==3.13.*"
# dependencies = ["dggal", "numpy>=1.24", "pyarrow>=15"]
# ///
"""Generate random DGGAL cell sets (ISEA7H, IVEA7H, ...) -> Parquet.

One file per system in dggal_common.DGGAL_SYSTEMS per level 0..target (coarse
levels with <= N cells are enumerated in full, finer ones sampled; schema +
write logic in _common.py). Reuses the parent dggal_common.py for the DGGAL
setup + registry, so adding a system is still one row there.

dggal ships an arch-broken macOS arm64 wheel, so on Apple Silicon this script
re-execs itself under an x86_64 (Rosetta) Python 3.13 (where the wheel is
self-consistent) — see the guard below. So it runs the same everywhere, with no
--python flag:

    uv run scripts/dggs/cells/gen_dggal.py

No CLI args (project convention) — edit the constants below in place.
"""

import os
import platform
import sys
from pathlib import Path

# On Apple Silicon, dggal's arm64 wheel bundles x86_64 dylibs and won't load, so
# re-exec under an x86_64 (Rosetta) Python 3.13 before importing it — keeping the
# command a plain `uv run` on every platform. Linux wheels are correct (no-op
# there). The env guard breaks the re-exec loop; under Rosetta machine() reports
# 'x86_64', so the condition is false on the second pass anyway.
if (sys.platform == 'darwin' and platform.machine() == 'arm64'
        and not os.environ.get('_DGGAL_ROSETTA')):
    os.environ['_DGGAL_ROSETTA'] = '1'
    os.execvp('uv', ['uv', 'run', '--python', 'cpython-3.13-macos-x86_64',
                     os.path.abspath(__file__)])

# dggal_common.py lives one level up (scripts/dggs/); reuse its DGGAL glue.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import dggal_common as dc  # noqa: E402

import _common  # noqa: E402

# Target resolution, N_BIG/N_SMALL, and SEED are pipeline config in _common.py.


if __name__ == '__main__':
    for name, row in dc.DGGAL_SYSTEMS.items():
        ad = dc.Adapter(row['cls'])
        # big: 0..target (_common.TARGET_RES); small: 0..the grid's finest level.
        _common.generate_big_small(
            name, _common.TARGET_RES[name], ad.max_level(),
            _common.N_BIG, _common.N_SMALL, _common.SEED,
            latlng_to_cell=lambda res, lat, lng, _ad=ad: _ad.zone_at(res, lng, lat),
            cid_str=ad.cid_str,
            cell_boundary=ad.ring_latlng,
            count_at=lambda res, _ad=ad: _ad.count(res),
            enumerate_at=lambda res, _ad=ad: _ad.enumerate(res),
        )

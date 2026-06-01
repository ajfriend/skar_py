# TODO

- [ ] **Backport the ~4s test loop to `sparea_py`.** sparea_py has the same old
  `reinstall: uv cache clean + uv sync --reinstall-package` + `test: reinstall
  ci-test` setup, so it pays the slow/occasionally-stalling rebuild. Port the
  final config from skar_py — note `no-build-isolation` must be **local only**
  (a global `[tool.uv]` setting breaks `uv build`/cibuildwheel in CI, which we
  hit):
  - `pyproject.toml`: a **non-default `build` group** with the backend
    (`meson-python`, `ninja`, `cython`, `ziglang`); keep `dev` to test deps.
    Do **not** add `[tool.uv] no-build-isolation-package`.
  - `justfile`: `reinstall: uv sync --reinstall-package sparea
    --no-build-isolation-package sparea --group build`; drop `uv cache clean`;
    `ci-test` uses `uv run --no-sync`. Stay non-editable.
  - `.github/workflows`: set `UV_NO_EDITABLE: "1"` on the test job; if the
    suite imports matplotlib, make it optional (skip plot tests when absent) so
    the minimal wheel/sdist test envs pass.
  - Also bump `mlugg/setup-zig@v1 -> @v2` (v1 404s on the moved Zig mirror).

- [ ] **Try a Rust-backed A5 binding in the DGGS survey.** A5 cell generation
  dominates the wall-clock of `scripts/dggs/survey.py` / `dggs_survey.ipynb`
  (~2.5s of ~2.7s for A5 at N=5000; H3 and S2 are ~0.1s each). The current
  `pya5` package looks pure-Python; if there's a Rust-backed A5 binding (the
  a5geo project has a Rust impl), swapping it in could make A5 as fast as the
  others and cut the survey time several-fold. Just a dependency swap in the
  `dggs`/`lab` groups if the API matches.

## Resolved

- **Fast non-editable test path** — done. `just test` is now ~4s, non-editable,
  no editable rebuild-on-import hook. Findings:
  - Zig was never the bottleneck (cold ReleaseFast ~4s, warm ~0.07s; Release ≈
    Debug). The rebuild *machinery* was.
  - uv does **not** rebuild a non-editable install on source change, so an
    explicit `--reinstall-package` per `just test` is genuinely required.
  - The ~2-minute hang was **not** the build: it was `uv cache clean` waiting
    on uv's global cache lock (`Timeout (300s) ... lock on ~/.cache/uv`) while
    another uv process held it. Dropping `uv cache clean` removes that failure
    mode and discards nothing useful.
  - `no-build-isolation` + build backend in the venv is what gets a rebuild to
    ~4s instead of ~9s (no re-staging of ziglang/meson/cython each time).
  - An editable install reaches ~0.5s but adds meson-python's rebuild-on-import
    hook (fragile `ninja` path under build isolation); not worth it for a ~4s
    gap.

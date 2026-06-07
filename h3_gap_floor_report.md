# skar convergence stall on a band of H3 resolutions (gap floors at ~1.7e-6)

**Status:** RESOLVED in skar_zig **v0.2.0** (2026-06-07) — the hypothesis below
was correct: the fix lowered `ACTIVE_THRESH` from 1e-6 to 1e-12. skar_py now
pins v0.2.0; the r7–r10 band converges at the strict `gap_tol = 1e-6` default
(verified: reproducer cell `899f4d0cd47ffff` now converges, AR 1.019514). The
diagnosis below is kept as the historical write-up.
**Found:** 2026-06-06, while moving `scripts/dggs/survey.py` from finest
resolution to H3 res 9 (a commonly-used resolution).
**Summary:** A non-trivial fraction (~1–3%) of random H3 cells in the
**r7–r10** band fail to certify the strict default `gap_tol = 1e-6`: the
duality gap stalls at ~1.7e-6 and never drops below 1e-6, no matter how many
outer iterations it runs. This looks **algorithmic, not an f64 precision
floor** — the opposite of the documented finest-resolution S2/A5 behaviour.
The aspect ratio itself is accurate regardless; only the certificate is stuck.

## Why this is surprising

`src/zig/.../tests/dggs_dnc_test.zig` already documents and pins the
finest-resolution DNCs for S2/A5: those are a genuine **f64 gap floor**
(sub-metre scatters at an O(1) point, κ(A) ~ 1e9, floor ~3.4e-4), and it
explicitly notes **"H3 r15 hits no f64 gap floor so it converges at the strict
1e-6 default"** — with a canary that an H3 r15 cell converges in **1 outer
iteration**.

So we expected coarser H3 (larger, *better*-conditioned cells) to be at least
as easy as r15. Instead a band of mid resolutions stalls.

## Evidence

### 1. It's a mid-resolution band, vanishing at the finest end

N = 10000 uniform-random cells per resolution, seed `0xC0FFEE`,
`gap_tol = 1e-6`, default `max_outer = 100`:

| H3 res | DNC rate | max outer iters (converged) | DNC gap floor (min/med/max) |
|--------|---------:|----------------------------:|-----------------------------|
| r15 (finest) | 0.0% | 1 | — |
| r12 | 0.0% | 1 | — |
| r10 | 3.0% | 7 | 1.27e-6 / 1.72e-6 / 2.18e-6 |
| r9  | 3.4% | 11 | 1.28e-6 / 1.70e-6 / 2.18e-6 |
| r8  | 1.3% | 4 | 1.27e-6 / 1.70e-6 / 2.15e-6 |
| r7  | 0.6% | 10 | 1.27e-6 / 1.75e-6 / 2.18e-6 |
| r5  | 0.0% | 10 | — |

A real f64 precision floor scales with κ(A) and would get **monotonically
worse** as cells shrink. Here it gets *better* toward the finest end (r12/r15
= 0%), peaks around r9–r10, and disappears at r5. That is not a precision
floor.

### 2. The stall value is scale-independent

Every affected resolution floors at the *same* ~1.27e-6–2.18e-6 band,
independent of cell size (which varies by ~50x across r7–r10). A κ-driven
floor would move with scale; this does not. The fixed ~1e-6-ish value sits
suspiciously next to a hard-coded solver constant (see Hypothesis).

### 3. Chasing 1e-6 makes the gap drift *up* (oscillation near the floor)

Single reproducer cell (see fixture below), same input, varying only
`gap_tol`:

| gap_tol | status | final gap | aspect ratio |
|---------|--------|-----------|--------------|
| 1e-6 | **did_not_converge** | **1.695e-6** | 1.01951364 (from sigma) |
| 2e-6 | converged | **1.341e-6** | 1.01951390 |
| 1e-5 | converged | 1.341e-6 | 1.01951390 |
| 1e-3 | converged | 1.341e-6 | 1.01951390 |

The cell can reach gap **1.341e-6**. But when asked for 1e-6 it runs the full
100 outer iterations and ends *worse*, at 1.695e-6 — the iterate wanders
upward instead of settling. So 1e-6 is unreachable for this cell while
~1.34e-6 is reachable; pushing past the achievable point degrades the
certificate. Raising `max_outer` to 2000 does **not** reduce the DNC count —
confirming it stalls/oscillates rather than slowly converging.

The aspect ratio is identical to ~7 significant figures across all tolerances
(1.0195139), so the *answer* callers want is correct; only the duality-gap
certificate is affected.

## Reproducer

H3 cell `899f4d0cd47ffff` (res 9), a near-circular hexagon (AR ≈ 1.0195).

Python:

```python
import h3, skar
v = skar.to_vec3(h3.cell_to_boundary('899f4d0cd47ffff'), geo='latlng')
skar.solve(v, geo='vec3').status                  # 'did_not_converge' (gap 1.695e-6)
skar.solve(v, geo='vec3', gap_tol=1e-5).status    # 'converged'  (gap 1.341e-6)
skar.solve(v, geo='vec3', gap_tol=1e-5).aspect_ratio   # 1.0195139048031325
```

Zig-ready fixture (unit vec3, same style as `dggs_dnc_test.zig`'s `A5_CELL`/
`S2_CELL`):

```zig
// H3 r9 cell 899f4d0cd47ffff — near-circular hexagon (AR ~1.0195). DNCs at
// the strict 1e-6 default (gap stalls ~1.7e-6); converges at >=2e-6.
const H3_R9_CELL = [_][3]f64{
    .{ -0.8586175701975843, 0.28761239723198995, -0.42432885490673883 },
    .{ -0.8586271933201559, 0.28762660191847433, -0.42429975342908594 },
    .{ -0.8586197375801148, 0.2876590246563569,  -0.42429286085392487 },
    .{ -0.8586026585975493, 0.2876772430738286,  -0.42431506980858175 },
    .{ -0.8585930353179254, 0.2876630384891841,  -0.42434417162336724 },
    .{ -0.8586004911779209, 0.28763061538522544, -0.42435106414636176 },
};
```

## Hypothesis (to test in skar_zig)

The fixed ~1e-6 stall lines up with `src/config.zig`:

```zig
pub const algo = struct {
    /// Certificate active-set cutoff: weights below this are dropped
    /// from `Info.cert`.
    pub const ACTIVE_THRESH: f64 = 1e-6;
    ...
};
```

If a binding constraint for these near-circular cells carries a dual weight
just **below** `ACTIVE_THRESH = 1e-6`, it gets dropped from the certificate,
leaving the computed gap stuck at ~1e-6 even though the primal/dual iterates
are essentially optimal. That would explain (a) the scale-independent value
and (b) why the gap can't be driven below ~1.3e-6.

A second, possibly related factor: `config.zig` flags near-isotropic `M`
(exactly the near-circular hex / DGGS case) as numerically delicate —
`PRECOND_COND_MIN` / `AXIS_WARMUP` exist specifically because the
preconditioner "adds sub-ULP direction noise that interacts badly with
damping after Newton polish" for these inputs. The oscillation in Evidence 3
(gap drifting 1.34e-6 -> 1.70e-6) is consistent with a damping/step-control
interaction once the iterate is already at the achievable optimum.

### Suggested investigation steps

1. Instrument `newton.zig` / the certificate assembly to dump, for cell
   `899f4d0cd47ffff`, the dual weights and the per-outer-iteration gap. Check
   whether a weight sits just under `ACTIVE_THRESH` and whether the gap
   oscillates after the first few iterations.
2. Try lowering `ACTIVE_THRESH` (e.g. 1e-9) and/or computing the certified gap
   *before* the active-set cutoff; see if the r7–r10 band converges at 1e-6
   without harming the finest-resolution S2/A5 guards.
3. Add an early-stop / best-iterate-tracking so chasing an unreachable
   tolerance returns the best achieved gap rather than drifting upward.
4. Extend `dggs_dnc_test.zig` with the H3 r9 fixture above as a regression
   once fixed (expect convergence at 1e-6), and keep an H3-band canary.

## Impact / current workaround

`scripts/dggs/survey.py` runs at the H3-r9-matched resolutions (H3 r9, S2 L15,
A5 r14) and uses `GAP_TOL = 1e-5` — comfortably above the ~2.2e-6 stall, so
all three systems converge with complete distributions, and still 100x
stricter than the old finest-resolution `1e-3`. The aspect ratios are accurate
regardless of the certificate, so this does not affect the survey's results.
The only reason to fix the solver is to honour the strict `1e-6` default on
these common, well-conditioned cells.

## Kickoff prompt (for a new session in the skar_zig repo)

> There's a convergence issue in the solver. Read this file
> (`h3_gap_floor_report.md`) first — it has the full evidence and a
> reproducer.
>
> **Symptom:** a band of H3 DGGS resolutions (r7–r10) has ~1–3% of cells whose
> duality-gap certificate stalls at ~1.7e-6 and never reaches the strict
> `gap_tol = 1e-6`. Reproducer is the `H3_R9_CELL` vec3 fixture above (H3 r9
> cell `899f4d0cd47ffff`, AR ≈ 1.0195): `.did_not_converge` at the 1e-6
> default, converges at ≥2e-6.
>
> **Why I think it's algorithmic, not an f64 floor** (verify, don't take on
> faith): it's a *mid*-resolution band — the finest cells (r12/r15) converge in
> 1 outer iteration; the stall value is *scale-independent* (~1.7e-6 across
> r7–r10 regardless of cell size); and chasing 1e-6 makes the gap drift
> *upward* (1.34e-6 → 1.70e-6) while raising `max_outer` to 2000 changes
> nothing. A κ-driven f64 floor would do the opposite.
>
> **Hypothesis to test:** `src/config.zig` has `algo.ACTIVE_THRESH = 1e-6` (the
> certificate active-set weight cutoff). If a binding constraint for these
> near-circular cells carries a dual weight just under 1e-6, it gets dropped
> from the certificate, flooring the computed gap near 1e-6.
>
> **What I want:**
> 1. Reproduce in Zig directly with the fixture. Instrument the per-outer-
>    iteration gap and the dual weights for this cell. Confirm or refute the
>    `ACTIVE_THRESH` hypothesis and identify the actual mechanism.
> 2. Only then propose a fix. Candidates: lower `ACTIVE_THRESH`, compute the
>    certified gap *before* the active-set cutoff, and/or track the best iterate
>    so an unreachable target returns the best achieved gap instead of drifting.
>
> **Guardrails:**
> - Do **not** "fix" this by loosening the default tolerance or making the gap
>   check less strict — the goal is to genuinely certify 1e-6 on these
>   well-conditioned cells.
> - `tests/dggs_dnc_test.zig` documents and pins the finest-resolution S2/A5
>   DNCs, which are a *correct* f64 floor. Your change must **not** make those
>   cells falsely "converge" at 1e-6. Run that test.
> - Add a regression test (extend `dggs_dnc_test.zig`) using the H3 r9 fixture:
>   it should converge at 1e-6 once fixed, plus an H3-band canary.

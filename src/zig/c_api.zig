//! C ABI shim for the skar Python bindings. Marshals between C
//! scalars / out-params and the upstream `skar.solve` API; all real
//! work happens in the `skar` Zig package.
//!
//! Minimal surface: it exposes the headline result of a solve —
//! status, aspect ratio, cone axis, eigenvalues, gap — but not the
//! variable-length active-set certificate (indices/lambdas). The cert
//! is freed here via `outcome.deinit()` before returning; surfacing it
//! across the C boundary would mean a caller-owned buffer protocol and
//! is deferred until a consumer needs it.

const std = @import("std");
const skar = @import("skar");

// Return codes for the `skar_solve` call itself: 0 = "ran, status is
// in out_status"; non-zero = "could not run". Mirrors the upstream
// errors-vs-outcome split (see src/api.zig in skar_zig).
pub const SKAR_OK: c_int = 0;
pub const SKAR_INSUFFICIENT_POINTS: c_int = 1;
pub const SKAR_INVALID_TOLERANCE: c_int = 2;
pub const SKAR_COPLANAR_INPUT: c_int = 3;
pub const SKAR_OUT_OF_MEMORY: c_int = 4;
pub const SKAR_INTERNAL: c_int = 5;
pub const SKAR_INVALID_METHOD: c_int = 6;

// `method` in-param values, mapping onto `skar.Method`. SKAR_METHOD_AUTO
// is upstream's alias for its recommended method (`Method.recommended`);
// `out_method` reports the concrete path that produced the outcome.
pub const SKAR_METHOD_ALTERNATING: c_int = 0;
pub const SKAR_METHOD_TRUST: c_int = 1;
pub const SKAR_METHOD_AUTO: c_int = 2;

// Values written to `out_status` on SKAR_OK — which Outcome variant
// the solver produced.
pub const SKAR_STATUS_CONVERGED: c_int = 0;
pub const SKAR_STATUS_INFEASIBLE: c_int = 1;
pub const SKAR_STATUS_DID_NOT_CONVERGE: c_int = 2;

/// The `out_method` value for a diagnostics union: which solver path
/// produced the outcome (the union tag).
fn pathTag(diag: skar.Diagnostics) c_int {
    return switch (diag) {
        .alternating => SKAR_METHOD_ALTERNATING,
        .trust => SKAR_METHOD_TRUST,
    };
}

/// C ABI: solve the spherical aspect-ratio problem for a point set.
///
/// `pts_buf` is an interleaved `(N, 3)` row-major buffer
/// `[x0, y0, z0, ...]` — what numpy hands us for a `(N, 3)` float64
/// array. `[3]f64` is exactly that layout, so we reinterpret the
/// pointer with no copy. `n_hull`, `gap_tol`, `coplanarity_tol`, and
/// `max_outer` map straight onto `skar.SolveOptions`.
///
/// On SKAR_OK, `out_status` selects which outputs are meaningful:
///   - converged:        sigma, q, gap, outer_iters
///   - did_not_converge: sigma, q, gap, outer_iters (uncertified)
///   - infeasible:       residual
/// Outputs not meaningful for the variant are left as NaN / 0. The
/// cone axis is column 0 of `q` (q[0], q[3], q[6]) and the aspect ratio
/// is sigma[2]/sigma[1], so neither is returned separately.
/// `out_outer_iters` is `Diagnostics.totalIters()` upstream — the
/// outer-iteration count on the alternating path; opening rounds +
/// trust-region iterations + re-certification attempts on the trust
/// path. `out_method` reports which path produced the outcome
/// (SKAR_METHOD_ALTERNATING/TRUST; under SKAR_METHOD_AUTO that is the
/// method the alias resolved to; -1 for infeasible, which carries no
/// path tag).
pub export fn skar_solve(
    pts_buf: [*]const f64,
    n: usize,
    gap_tol: f64,
    n_hull: c_int,
    coplanarity_tol: f64,
    max_outer: c_uint,
    method: c_int,
    out_status: *c_int,
    out_sigma: *[3]f64,
    out_q: *[9]f64,
    out_gap: *f64,
    out_outer_iters: *c_uint,
    out_residual: *f64,
    out_method: *c_int,
) c_int {
    const X: []const [3]f64 = @as([*]const [3]f64, @ptrCast(pts_buf))[0..n];

    const opts: skar.SolveOptions = .{
        .gap_tol = gap_tol,
        .n_hull = @intCast(n_hull),
        .coplanarity_tol = coplanarity_tol,
        .max_outer = @intCast(max_outer),
        .method = switch (method) {
            SKAR_METHOD_ALTERNATING => .alternating,
            SKAR_METHOD_TRUST => .trust,
            SKAR_METHOD_AUTO => .auto,
            else => return SKAR_INVALID_METHOD,
        },
    };

    const nan = std.math.nan(f64);
    out_sigma.* = .{ nan, nan, nan };
    // Q is row-major (out_q[r*3 + c] = Q(r, c)); column i is the unit
    // eigenvector paired with sigma[i], and column 0 is the cone axis.
    out_q.* = .{nan} ** 9;
    out_gap.* = nan;
    out_outer_iters.* = 0;
    out_residual.* = nan;
    out_method.* = -1;

    var outcome = skar.solve(std.heap.c_allocator, X, opts) catch |err| switch (err) {
        error.InsufficientPoints => return SKAR_INSUFFICIENT_POINTS,
        error.InvalidTolerance => return SKAR_INVALID_TOLERANCE,
        error.CoplanarInput => return SKAR_COPLANAR_INPUT,
        error.OutOfMemory => return SKAR_OUT_OF_MEMORY,
        // SolveError variants (NegativeDualityGap / NegativeEigenvalue /
        // SingularMoment): internal-correctness bugs in the library.
        else => return SKAR_INTERNAL,
    };
    defer outcome.deinit();

    switch (outcome) {
        .converged => |c| {
            out_status.* = SKAR_STATUS_CONVERGED;
            out_sigma.* = c.sigma;
            out_q.* = c.Q.m;
            out_gap.* = c.gap;
            out_outer_iters.* = c.diag.totalIters();
            out_method.* = pathTag(c.diag);
        },
        .infeasible => |inf| {
            out_status.* = SKAR_STATUS_INFEASIBLE;
            out_residual.* = inf.residual;
        },
        .did_not_converge => |d| {
            out_status.* = SKAR_STATUS_DID_NOT_CONVERGE;
            out_sigma.* = d.sigma;
            out_q.* = d.Q.m;
            out_gap.* = d.gap;
            out_outer_iters.* = d.diag.totalIters();
            out_method.* = pathTag(d.diag);
        },
    }
    return SKAR_OK;
}

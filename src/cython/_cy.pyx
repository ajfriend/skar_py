# cython: language_level=3
"""Cython binding for skar — internal; called from `skar.solve`.

Compiled by meson (driven by meson-python); links against the Zig
static archive libskar.{a,lib}. Exposed as `skar._cy`.
"""

cdef extern from *:
    """
    int skar_solve(const double *pts, size_t n,
                   double gap_tol, int n_hull, double coplanarity_tol,
                   unsigned int max_outer, int method,
                   int *out_status,
                   double *out_sigma, double *out_q, double *out_gap,
                   unsigned int *out_outer_iters, double *out_residual,
                   int *out_method);
    """
    int skar_solve(const double *pts, size_t n,
                   double gap_tol, int n_hull, double coplanarity_tol,
                   unsigned int max_outer, int method,
                   int *out_status,
                   double *out_sigma, double *out_q, double *out_gap,
                   unsigned int *out_outer_iters, double *out_residual,
                   int *out_method)


# Status strings indexed by the SKAR_STATUS_* code that c_api.zig
# writes to out_status. Owned here, at the C boundary, so the raw
# integer codes never leak into the Python wrapper.
_STATUS = ('converged', 'infeasible', 'did_not_converge')

# Solver-path strings indexed by the SKAR_METHOD_* code; also the
# in-param encoding ('auto' = index 2 is input-only — out_method always
# reports a concrete path, or -1 (None) for infeasible).
_METHOD = ('alternating', 'trust', 'auto')


def solve(double[:, ::1] pts not None, double gap_tol, int n_hull,
          double coplanarity_tol, unsigned int max_outer, int method):
    if pts.shape[1] != 3:
        raise ValueError('pts must be a 2-D array of shape (N, 3)')

    cdef int status
    cdef double gap, residual
    cdef double sigma[3]
    cdef double q[9]
    cdef unsigned int outer_iters
    cdef int out_method
    cdef int err = skar_solve(
        &pts[0, 0], pts.shape[0], gap_tol, n_hull, coplanarity_tol, max_outer,
        method,
        &status, &sigma[0], &q[0], &gap, &outer_iters, &residual, &out_method,
    )

    if err == 1:
        raise ValueError('skar: need at least 3 points to define a cone')
    if err == 2:
        raise ValueError('skar: tolerances must be finite and positive')
    if err == 3:
        raise ValueError(
            'skar: input is near-coplanar (points lie ~on a great circle); '
            'pass coplanarity_tol<=0 to bypass this check'
        )
    if err == 4:
        raise MemoryError('skar: out of memory')
    if err == 5:
        raise RuntimeError(
            'skar: internal solver error (a PSD/duality invariant was '
            'violated beyond float noise) — please report it'
        )
    if err == 6:
        raise ValueError("skar: method must be 'alternating', 'trust', or 'auto'")
    if err != 0:
        raise RuntimeError(f'skar: unknown error code {err}')

    return (
        _STATUS[status],
        (sigma[0], sigma[1], sigma[2]),
        (q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7], q[8]),
        gap,
        outer_iters,
        residual,
        None if out_method < 0 else _METHOD[out_method],
    )

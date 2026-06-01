import math

import numpy as np
import pytest

import skar


# Octant: north pole + two equator points 90° apart. By 3-fold
# symmetry around (1,1,1)/√3 the tightest enclosing cone is circular,
# so the aspect ratio is 1.
OCTANT_DEG = np.array([
    [0.0,  0.0],   # (1, 0, 0)
    [0.0, 90.0],   # (0, 1, 0)
    [90.0, 0.0],   # (0, 0, 1)
])

OCTANT_RAD = np.radians(OCTANT_DEG)

OCTANT_XYZ = np.array([
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
])


@pytest.mark.parametrize('points,geo', [
    (OCTANT_DEG, 'latlng'),       # default = degrees
    (OCTANT_DEG, 'latlng_deg'),   # explicit degrees
    (OCTANT_RAD, 'latlng_rad'),   # explicit radians
    (OCTANT_XYZ, 'vec3'),
])
def test_octant_is_circular_cone(points, geo):
    r = skar.solve(points, geo=geo)
    assert r.status == 'converged'
    assert math.isclose(r.aspect_ratio, 1.0, abs_tol=1e-6)
    # Cone axis is (1,1,1)/√3.
    axis = np.array(r.axis)
    assert math.isclose(abs(axis @ (np.ones(3) / np.sqrt(3))), 1.0, abs_tol=1e-6)


def test_latlng_matches_vec3():
    a = skar.solve(OCTANT_DEG).aspect_ratio
    b = skar.solve(OCTANT_XYZ, geo='vec3').aspect_ratio
    assert math.isclose(a, b, rel_tol=1e-9)


def test_elongated_scatter_has_aspect_above_one():
    # A short, wide arc near the north pole: clearly non-circular, so
    # aspect ratio should be comfortably > 1.
    pts = np.array([
        [80.0, -40.0],
        [82.0,   0.0],
        [80.0,  40.0],
        [85.0,   0.0],
    ])
    r = skar.solve(pts)
    assert r.status == 'converged'
    assert r.aspect_ratio > 1.0


def test_primary_usage_is_a_plain_list_of_latlng_points():
    # The intended ergonomics: a list of (lat, lng) tuples, no numpy at
    # the call site.
    pts = [
        (0.0,  0.0),
        (0.0, 90.0),
        (90.0, 0.0),
    ]
    r = skar.solve(pts)
    assert r.status == 'converged'
    assert math.isclose(r.aspect_ratio, 1.0, abs_tol=1e-6)


def test_accepts_python_list_of_vec3():
    pts = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    r = skar.solve(pts, geo='vec3')
    assert r.status == 'converged'
    assert math.isclose(r.aspect_ratio, 1.0, abs_tol=1e-6)


def test_infeasible_when_no_hemisphere_contains_points():
    # Regular tetrahedron vertices span the sphere — no hemisphere
    # holds all four, so the problem is infeasible.
    pts = np.array([
        [1.0,  1.0,  1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ])
    pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    r = skar.solve(pts, geo='vec3')
    assert r.status == 'infeasible'
    assert r.aspect_ratio is None
    assert r.residual is not None
    assert r.residual < 1e-6


def test_invalid_geo():
    with pytest.raises(ValueError, match='geo must be'):
        skar.solve(OCTANT_DEG, geo='xyz')


def test_shape_mismatch():
    with pytest.raises(ValueError, match=r'\(N, 3\)'):
        skar.solve(OCTANT_DEG, geo='vec3')  # 2 cols for a 3-col geo


def test_shape_not_2d():
    with pytest.raises(ValueError):
        skar.solve(np.zeros(6))


def test_too_few_points():
    with pytest.raises(ValueError, match='at least 3'):
        skar.solve(OCTANT_XYZ[:2], geo='vec3')


# A short arc on the equator (a single great circle), contained well
# within a hemisphere: the tangent-plane projection is near-collinear,
# which is the CoplanarInput case.
COPLANAR_ARC_DEG = np.array([
    [0.0, -10.0],
    [0.0,  -3.0],
    [0.0,   4.0],
    [0.0,  10.0],
])


def test_coplanar_input_rejected():
    with pytest.raises(ValueError, match='coplanar'):
        skar.solve(COPLANAR_ARC_DEG)


def test_coplanar_check_can_be_bypassed():
    # Same near-collinear scatter, but bypass the check — the solver may
    # converge or not, but it must not raise CoplanarInput.
    r = skar.solve(COPLANAR_ARC_DEG, coplanarity_tol=0.0)
    assert r.status in {'converged', 'infeasible', 'did_not_converge'}


def test_invalid_tolerance():
    with pytest.raises(ValueError, match='tolerance'):
        skar.solve(OCTANT_XYZ, geo='vec3', gap_tol=float('nan'))

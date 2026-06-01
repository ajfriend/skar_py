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
    assert isinstance(r, skar.Converged)
    assert math.isclose(r.aspect_ratio, 1.0, abs_tol=1e-6)
    # Cone axis (first column of Q) is (1,1,1)/√3.
    axis = r.Q[:, 0]
    assert math.isclose(abs(axis @ (np.ones(3) / np.sqrt(3))), 1.0, abs_tol=1e-6)


def test_converged_exposes_orthonormal_Q():
    r = skar.solve(OCTANT_XYZ, geo='vec3')
    assert isinstance(r, skar.Converged)
    assert isinstance(r, skar.Outcome)  # native | union supports isinstance
    assert isinstance(r.Q, np.ndarray) and r.Q.shape == (3, 3)
    assert isinstance(r.sigma, np.ndarray) and r.sigma.shape == (3,)
    # Orthonormal and right-handed.
    assert np.allclose(r.Q.T @ r.Q, np.eye(3), atol=1e-9)
    assert math.isclose(np.linalg.det(r.Q), 1.0, abs_tol=1e-9)
    # sigma[0] is the structural axial eigenvalue 1/sqrt(3).
    assert math.isclose(r.sigma[0], 1.0 / math.sqrt(3.0), abs_tol=1e-9)
    # Aspect ratio is sigma[2] / sigma[1].
    assert math.isclose(r.aspect_ratio, r.sigma[2] / r.sigma[1], rel_tol=1e-12)


def test_Q_and_sigma_reconstruct_a_feasible_cone():
    # A = Q diag(sigma) Qᵀ; every input point must satisfy
    # ‖A·x‖ − b·x <= 0 (it sits inside the certified cone). The axis b
    # is simply Q[:, 0].
    pts = OCTANT_XYZ
    r = skar.solve(pts, geo='vec3')
    A = r.Q @ np.diag(r.sigma) @ r.Q.T
    b = r.Q[:, 0]
    viol = np.array([np.linalg.norm(A @ x) - b @ x for x in pts])
    assert np.all(viol <= 1e-6)


def test_to_vec3_converts_latlng_to_unit_vectors():
    # Octant (lat, lng) degrees → the identity basis on the sphere.
    X = skar.to_vec3([(0, 0), (0, 90), (90, 0)])
    assert isinstance(X, np.ndarray) and X.shape == (3, 3)
    assert np.allclose(X, np.eye(3), atol=1e-12)
    assert np.allclose(np.linalg.norm(X, axis=1), 1.0)


def test_to_vec3_vec3_is_passthrough():
    pts = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    assert np.allclose(skar.to_vec3(pts, geo='vec3'), pts)


def test_to_vec3_reads_geo_interface_as_lng_lat():
    geom = _Geo({'type': 'MultiPoint', 'coordinates': [[0.0, 0.0], [90.0, 0.0]]})
    # GeoJSON (lng, lat): (0,0)->(1,0,0); (90,0)->(0,1,0).
    assert np.allclose(skar.to_vec3(geom), [[1, 0, 0], [0, 1, 0]], atol=1e-12)


def test_solve_matches_explicit_to_vec3():
    pts = [(0, 0), (0, 90), (90, 0)]
    a = skar.solve(pts)
    b = skar.solve(skar.to_vec3(pts), geo='vec3')
    assert a.aspect_ratio == b.aspect_ratio


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
    assert isinstance(r, skar.Converged)
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
    assert isinstance(r, skar.Converged)
    assert math.isclose(r.aspect_ratio, 1.0, abs_tol=1e-6)


def test_accepts_python_list_of_vec3():
    pts = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    r = skar.solve(pts, geo='vec3')
    assert isinstance(r, skar.Converged)
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
    assert isinstance(r, skar.Infeasible)
    assert r.status == 'infeasible'
    assert r.residual < 1e-6
    # The type carries only its relevant field — there is no silent
    # None to misuse; converged-only fields simply don't exist.
    assert not hasattr(r, 'aspect_ratio')
    assert not hasattr(r, 'sigma')
    assert not hasattr(r, 'Q')


def test_did_not_converge_exposes_uncertified_diagnostics():
    # An unreachable gap tolerance (far below the f64 conditioning
    # floor) forces the solver to exhaust its iteration budget.
    pts = np.array([
        [80.0, -40.0],
        [82.0,   0.0],
        [80.0,  40.0],
        [85.0,   0.0],
    ])
    r = skar.solve(pts, gap_tol=1e-300, max_outer=5)
    assert isinstance(r, skar.DidNotConverge)
    assert r.status == 'did_not_converge'
    assert r.Q.shape == (3, 3) and r.sigma.shape == (3,)
    assert r.outer_iters == 5
    # Uncertified: the aspect_ratio accessor is deliberately withheld.
    assert not hasattr(r, 'aspect_ratio')


class _Geo:
    """Minimal stand-in for shapely/geojson/h3 — anything exposing the
    `__geo_interface__` protocol. Avoids a real geometry dependency."""

    def __init__(self, mapping):
        self._mapping = mapping

    @property
    def __geo_interface__(self):
        return self._mapping


def test_geo_interface_polygon():
    # A "square" of points around the north pole, as a closed GeoJSON
    # (lng, lat) ring. By 4-fold symmetry the enclosing cone is
    # circular, so aspect ~ 1.
    poly = _Geo({
        'type': 'Polygon',
        'coordinates': [[
            [0.0, 80.0], [90.0, 80.0], [180.0, 80.0], [270.0, 80.0],
            [0.0, 80.0],   # closing duplicate vertex
        ]],
    })
    r = skar.solve(poly)
    assert r.status == 'converged'
    assert math.isclose(r.aspect_ratio, 1.0, abs_tol=1e-3)


def test_geo_interface_swaps_lng_lat():
    # GeoJSON is (lng, lat); confirm the geo-interface path matches
    # feeding the same points manually swapped to (lat, lng).
    positions = [[10.0, 80.0], [-20.0, 82.0], [40.0, 85.0], [5.0, 78.0]]
    via_geo = skar.solve(_Geo({'type': 'MultiPoint', 'coordinates': positions}))
    manual = skar.solve(
        [(lat, lng) for lng, lat in positions], geo='latlng_deg'
    )
    assert via_geo.status == manual.status == 'converged'
    assert math.isclose(via_geo.aspect_ratio, manual.aspect_ratio, rel_tol=1e-9)


def test_geo_interface_ignores_geo_argument():
    # geo='vec3' would be nonsense for (lng, lat) degrees, but the
    # geometry defines its own convention so `geo` is ignored.
    positions = [[10.0, 80.0], [-20.0, 82.0], [40.0, 85.0], [5.0, 78.0]]
    geom = _Geo({'type': 'LineString', 'coordinates': positions})
    assert skar.solve(geom, geo='vec3').status == 'converged'


def test_geo_interface_multipolygon():
    geom = _Geo({
        'type': 'MultiPolygon',
        'coordinates': [
            [[[0.0, 80.0], [90.0, 80.0], [180.0, 80.0], [0.0, 80.0]]],
            [[[270.0, 82.0], [300.0, 84.0], [330.0, 83.0], [270.0, 82.0]]],
        ],
    })
    assert skar.solve(geom).status == 'converged'


def test_geo_interface_feature_unwrap():
    feature = _Geo({
        'type': 'Feature',
        'properties': {},
        'geometry': {
            'type': 'MultiPoint',
            'coordinates': [[10.0, 80.0], [-20.0, 82.0], [40.0, 85.0]],
        },
    })
    assert skar.solve(feature).status == 'converged'


def test_geo_interface_unsupported_type():
    pt = _Geo({'type': 'Point', 'coordinates': [10.0, 80.0]})
    with pytest.raises(ValueError, match='unsupported __geo_interface__'):
        skar.solve(pt)


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

"""Input normalization: turn arbitrary point input into the `(N, 3)`
unit-vector array the solver consumes.
"""

from typing import Literal, get_args

import numpy as np

Geo = Literal['latlng', 'latlng_deg', 'latlng_rad', 'vec3']
_GEO = frozenset(get_args(Geo))


def _ring(ring):
    # A GeoJSON polygon ring is closed — the last position repeats the
    # first (RFC 7946 §3.1.6). Drop it so that vertex isn't counted twice.
    return list(ring[:-1])


def _geo_interface_positions(gi):
    """Flatten a GeoJSON-like mapping (the value of an object's
    ``__geo_interface__``) into a flat list of positions. GeoJSON
    positions are ``[lng, lat, ...]``; the caller swaps to skar's
    ``(lat, lng)``. Polygons contribute their exterior ring only."""
    match gi.get('type'):
        case 'Feature':
            return _geo_interface_positions(gi['geometry'])
        case 'MultiPoint' | 'LineString':
            return list(gi['coordinates'])
        case 'Polygon':
            return _ring(gi['coordinates'][0])
        case 'MultiPolygon':
            return [pos for poly in gi['coordinates'] for pos in _ring(poly[0])]
        case other:
            raise ValueError(
                f'skar: unsupported __geo_interface__ geometry type {other!r}; '
                f'expected MultiPoint, LineString, Polygon, or MultiPolygon'
            )


def to_vec3(points, *, geo: Geo = 'latlng') -> np.ndarray:
    """Convert `points` to the `(N, 3)` unit-vector array the solver
    consumes.

    This is exactly the conversion `solve` runs on its input; call it
    directly to see how your points land on the unit sphere.

    Args:
        points: a sequence of points — typically a list or tuple of
            `(lat, lng)` pairs (or `(x, y, z)` triples for
            `geo='vec3'`); each element is one point. A NumPy array is
            also accepted, interpreted as a 2-D `(N, k)` array whose
            **rows are points** and columns are coordinates (`k` = 2
            for the `latlng` family, 3 for `'vec3'`).

            An object exposing `__geo_interface__` (shapely, geojson,
            h3 `LatLngPoly`/`LatLngMultiPoly`, …) is also accepted: its
            vertices are extracted and read as GeoJSON `(lng, lat)`
            degrees. Supported geometry types are `MultiPoint`,
            `LineString`, `Polygon` (exterior ring), `MultiPolygon`,
            and a `Feature` wrapping one of those. For such inputs the
            `geo` argument is **ignored** — the geometry defines its own
            convention.
        geo: input convention.
            `'latlng'` (default) and `'latlng_deg'`: each point is
                `(lat, lng)` in **degrees** — matching h3's convention.
            `'latlng_rad'`: each point is `(lat, lng)` in radians.
            `'vec3'`: each point is a unit `(x, y, z)` on the sphere.

    Returns:
        A C-contiguous `(N, 3)` float64 array of points on the unit
        sphere. `latlng` inputs are converted to `(x, y, z)`; `'vec3'`
        rows are returned as-is (assumed already unit length).

    Raises:
        ValueError: invalid `geo`, mismatched shape, or unsupported
            `__geo_interface__` geometry.
    """
    # An object exposing __geo_interface__ defines its own convention:
    # GeoJSON (lng, lat) degrees. Extract its vertices, swap to skar's
    # (lat, lng), and ignore `geo`.
    if hasattr(points, '__geo_interface__'):
        positions = _geo_interface_positions(points.__geo_interface__)
        points = [(p[1], p[0]) for p in positions]
        geo = 'latlng_deg'

    if geo not in _GEO:
        raise ValueError(f'geo must be one of {sorted(_GEO)}, got {geo!r}')

    arr = np.ascontiguousarray(points, dtype=np.float64)
    cols = 3 if geo == 'vec3' else 2
    if arr.ndim != 2 or arr.shape[1] != cols:
        raise ValueError(
            f'points must be a 2-D array with shape (N, {cols}) for '
            f'geo={geo!r}, got shape {arr.shape}'
        )

    if geo == 'vec3':
        return arr

    lat, lng = arr[:, 0], arr[:, 1]
    if geo != 'latlng_rad':
        lat, lng = np.radians(lat), np.radians(lng)
    cl = np.cos(lat)
    # column_stack returns a fresh C-contiguous (N, 3) f64.
    return np.column_stack([cl * np.cos(lng), cl * np.sin(lng), np.sin(lat)])

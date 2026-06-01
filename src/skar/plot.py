"""Plotting helpers for an enclosing cone.

`project_to_cone` (pure NumPy) gnomonic-projects points at the cone axis and
returns the cross-section ellipse — the math behind the plots. `plot_cone`
draws a geometry's outline with that ellipse and needs matplotlib (install it
directly or via the `plot` extra, `pip install skar[plot]`).
"""

import numpy as np

from .convert import _geo_interface_rings, to_vec3
from .outcomes import Converged

EARTH_RADIUS_M = 6_371_008.8  # mean Earth radius; the default projection scale


def project_to_cone(result, points, *, up=(0.0, 0.0, 1.0), scale=EARTH_RADIUS_M):
    """Gnomonic-project unit-vector `points` onto the tangent plane at the cone
    axis, and return the cross-section ellipse's semi-axes — both in the same
    scaled frame, with the major axis horizontal.

    Working in skar's eigenbasis ``Q`` (axis ``Q[:, 0]``, tangent eigenvectors
    ``Q[:, 1:]``), the cross-section is axis-aligned with semi-axes
    ``sqrt(2/3) / sigma[1:]``, so no extra rotation is needed.

    Args:
        result: a `Converged` outcome from `solve`.
        points: an ``(N, 3)`` array of unit vectors on the sphere (e.g. from
            `to_vec3`).
        up: reference direction (a 3-vector in the cone's space) placed in the
            upper half via a 180° flip — which keeps the major axis horizontal.
            Defaults to the north pole `(0, 0, 1)`; pass `None` to skip.
        scale: multiply projected coordinates. Defaults to the Earth's radius
            in metres.

    Returns:
        ``(xy, semi)`` — `xy` is the ``(N, 2)`` projected coordinates; `semi` is
        the ``(major, minor)`` ellipse semi-axis lengths. The enclosing-cone
        cross-section is ``(xy[:, 0] / semi[0])**2 + (xy[:, 1] / semi[1])**2 = 1``,
        and ``semi[0] / semi[1]`` is the aspect ratio.

    Raises:
        TypeError: `result` is not a `Converged`.
    """
    if not isinstance(result, Converged):
        raise TypeError(
            f'project_to_cone needs a Converged result, got {type(result).__name__}'
        )

    b, U = result.Q[:, 0], result.Q[:, 1:]
    semi = np.sqrt(1.0 - result.sigma[0] ** 2) / result.sigma[1:] * scale

    flip = 1.0
    if up is not None:
        u = np.asarray(up, dtype=float)
        if ((u - (u @ b) * b) @ U)[1] < 0.0:  # up's tangent-plane y-component
            flip = -1.0

    points = np.asarray(points, dtype=float)
    xy = (points @ U) / (points @ b)[:, None] * flip * scale
    return xy, semi


def plot_cone(result, geometry, *, title=None, up=(0.0, 0.0, 1.0),
              scale=EARTH_RADIUS_M, size=7.0, ax=None):
    """Draw a geometry's outline gnomonic-projected at its enclosing-cone axis,
    with the cone's cross-section ellipse overlaid (major axis horizontal).

    Sensible defaults for a point set on the globe: the projection is scaled to
    metres, the panel is oriented north-up, and the axes are labelled — so a
    bare ``plot_cone(result, geometry)`` already produces a finished figure.

    Args:
        result: a `Converged` outcome from `solve` (the same geometry).
        geometry: an object exposing `__geo_interface__` (shapely / geopandas /
            geojson). Its rings (exterior + holes) are drawn separately.
        title: optional Axes title.
        up: orientation reference; see `project_to_cone`. Defaults to north-up.
        scale: projected-coordinate scale; see `project_to_cone`.
        size: side length in inches of the square figure created when `ax` is
            omitted.
        ax: matplotlib Axes to draw on; a new one is created if omitted.

    Returns:
        The matplotlib Axes.

    Raises:
        TypeError: `result` is not a `Converged`.
    """
    import matplotlib.pyplot as plt

    if not isinstance(result, Converged):  # fail before creating a figure
        raise TypeError(
            f'plot_cone needs a Converged result, got {type(result).__name__}'
        )

    if ax is None:
        _, ax = plt.subplots(figsize=(size, size))

    semi = None
    for ring in _geo_interface_rings(geometry.__geo_interface__):
        v = to_vec3(np.asarray(ring)[:, :2], geo='lonlat')
        xy, semi = project_to_cone(result, v, up=up, scale=scale)
        ax.plot(xy[:, 0], xy[:, 1], color='C0', lw=1.0)

    t = np.linspace(0.0, 2.0 * np.pi, 400)
    ax.plot(semi[0] * np.cos(t), semi[1] * np.sin(t), color='0.25', lw=1.5)
    ax.set_aspect('equal')
    ax.set_xlabel('major axis (m)')
    ax.set_ylabel('minor axis (m)')
    if title is not None:
        ax.set_title(title)
    return ax

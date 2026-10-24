"""Geometric functions of the section with exact gradients w.r.t. the design
variables.  Each function returns ``(values, jac)`` with ``values`` of shape
(m,) and ``jac`` of shape (m, ndv).

Linear functions (thickness, camber, box) use the parametrization's surface
operators; the polygon integrals (area, inertia) use the complex-step
derivative through the affine node map, exact to machine precision.
"""
import numpy as np


def _complex_step(fun, param, dv, h=1e-30):
    dv = np.asarray(dv, float)
    f0 = np.atleast_1d(fun(*param.coordinates(dv)).real)
    jac = np.zeros((len(f0), len(dv)))
    for k in range(len(dv)):
        d = dv.astype(complex)
        d[k] += 1j * h
        jac[:, k] = np.atleast_1d(fun(*param.coordinates(d))).imag / h
    return f0, jac


def section_properties(x, y):
    """Area, centroid and centroidal second moments of the closed polygon
    (counterclockwise nodes, the TE gap is closed by a straight segment).

    Returns dict: area, xc, yc, Ixx (about the centroidal axis parallel to the
    chord: flapwise bending), Iyy (about the centroidal axis normal to the
    chord), Ixy, and Ip = Ixx + Iyy."""
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    c = x * y1 - x1 * y
    A = 0.5 * np.sum(c)
    Sx = np.sum((y + y1) * c) / 6.0      # int y dA
    Sy = np.sum((x + x1) * c) / 6.0      # int x dA
    Ixx0 = np.sum((y * y + y * y1 + y1 * y1) * c) / 12.0
    Iyy0 = np.sum((x * x + x * x1 + x1 * x1) * c) / 12.0
    Ixy0 = np.sum((x * y1 + 2 * x * y + 2 * x1 * y1 + x1 * y) * c) / 24.0
    xc, yc = Sy / A, Sx / A
    Ixx = Ixx0 - A * yc * yc
    Iyy = Iyy0 - A * xc * xc
    Ixy = Ixy0 - A * xc * yc
    return dict(area=A, xc=xc, yc=yc, Ixx=Ixx, Iyy=Iyy, Ixy=Ixy, Ip=Ixx + Iyy)


def area(param, dv):
    return _complex_step(lambda x, y: section_properties(x, y)["area"], param, dv)


def inertia(param, dv, axis="xx"):
    """Centroidal second moment of area: axis in {'xx', 'yy', 'xy', 'p'}."""
    key = {"xx": "Ixx", "yy": "Iyy", "xy": "Ixy", "p": "Ip"}[axis]
    return _complex_step(lambda x, y: section_properties(x, y)[key], param, dv)


def section_modulus(param, dv):
    """Ixx / (max distance from the centroidal axis): bending stress criterion."""
    def fun(x, y):
        p = section_properties(x, y)
        d = np.sqrt((y - p["yc"]) ** 2 + 1e-12)       # smooth |.| (complex-step safe)
        # smooth max (KS) to keep it differentiable
        rho = 400.0
        dm = d.real.max()
        ks = dm + np.log(np.sum(np.exp(rho * (d - dm)))) / rho
        return p["Ixx"] / ks
    return _complex_step(fun, param, dv)


def thickness(param, dv, xs):
    c, R = param.thickness_op(np.atleast_1d(xs))
    return c + R @ dv, R


def camber(param, dv, xs):
    c, R = param.camber_op(np.atleast_1d(xs))
    return c + R @ dv, R


def max_thickness(param, dv, x_range=(0.02, 0.98), n=200, rho=2000.0):
    """Smooth (Kreisselmeier-Steinhauser) maximum of the thickness t(x)."""
    xs = np.linspace(*x_range, n)
    t, R = thickness(param, dv, xs)
    tm = t.max()
    w = np.exp(rho * (t - tm))
    ks = tm + np.log(w.sum()) / rho
    return np.array([ks]), ((w / w.sum()) @ R)[None, :]


def box(param, dv, x_range, y_range, n=15):
    """The section must enclose the rectangle x_range x y_range:
    returns g >= 0 with g = [y_upper(x_j) - y_top, y_bottom - y_lower(x_j)]."""
    xs = np.linspace(x_range[0], x_range[1], n)
    cu, Ru = param.upper_op(xs)
    cl, Rl = param.lower_op(xs)
    g = np.concatenate([cu + Ru @ dv - y_range[1], y_range[0] - (cl + Rl @ dv)])
    J = np.vstack([Ru, -Rl])
    return g, J


def le_radius(param, dv):
    """Leading-edge radius (B-spline parametrization: r is proportional to a^2)."""
    if not hasattr(param, "le_radius"):
        raise ValueError("le_radius is only available for BSplineAirfoil")
    r = param.le_radius(dv)
    J = np.zeros((1, len(dv)))
    J[0, 0] = 2.0 * r / dv[0]
    return np.array([r]), J


def te_angle(param, dv, dx=0.02):
    """Wedge angle (rad) between upper and lower surfaces over the last dx of chord."""
    xs = np.array([1.0 - dx, 1.0])
    t, R = thickness(param, dv, xs)
    val = (t[1] - t[0]) / dx            # negative: thickness decreasing to the TE
    ang = -2.0 * np.arctan(0.5 * val)
    dang = -2.0 / (1.0 + (0.5 * val) ** 2) * 0.5 / dx * (R[1] - R[0])
    return np.array([ang]), dang[None, :]


def smoothness(param, dv, side="both"):
    """Second differences of the B-spline control-point displacements from the
    initial shape (one value per interior control point and side).  Bounding
    them (lower=-s, upper=s) prevents wiggles when there are many control points."""
    if not hasattr(param, "smoothness_op"):
        raise ValueError("smoothness is only available for BSplineAirfoil")
    c, R = param.smoothness_op(side)
    return c + R @ dv, R


REGISTRY = {
    "smoothness": smoothness,
    "area": area,
    "inertia": inertia,
    "section_modulus": section_modulus,
    "thickness": thickness,
    "camber": camber,
    "max_thickness": max_thickness,
    "box": box,
    "le_radius": le_radius,
    "te_angle": te_angle,
}

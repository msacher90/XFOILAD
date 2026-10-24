"""Airfoil shape parametrizations, linear in their design variables.

Every parametrization maps a design vector ``dv`` to the XFOIL panel nodes
(counterclockwise: upper TE -> LE -> lower TE).  The map is affine:

    y_nodes = y0 + Jy @ dv        x_nodes = x0 + Jx @ dv

so the Jacobian is constant and the chain rule with the XFOIL adjoint
(gradient w.r.t. node coordinates) is exact and costless.

* ``BSplineAirfoil``: two clamped B-spline curves (upper / lower) joined at
  the leading edge with a vertical tangent; the y of the control points are
  the design variables, the x are fixed.  The first control point above /
  below the LE is shared (``a`` and ``-a``) which gives a curvature-continuous
  leading edge whose radius is driven by ``a``.
* ``CSTAirfoil``: Kulfan class/shape transformation, Bernstein coefficients
  as design variables.

Both take ``symmetric=True`` for a symmetric section: the lower surface is
the mirror of the upper one and only the upper-surface variables are design
variables (half as many).

``le_node``: with True (default of a non-symmetric section) the leading edge
is a panel node; with False (default of a symmetric section) it is replaced
by two nodes, one on each side, half way along the first panel.  At alpha = 0
the stagnation point of a symmetric section then lies in the middle of a
panel instead of on a node, where XFOIL's BL similarity station is singular.
"""
import numpy as np
from scipy.interpolate import BSpline


# --------------------------------------------------------------------------
# airfoil coordinates helpers
# --------------------------------------------------------------------------
def naca4(code="0012", n=200):
    """NACA 4-digit airfoil (closed form, finite TE), counterclockwise from upper TE."""
    m, p, t = int(code[0]) / 100.0, int(code[1]) / 10.0, int(code[2:]) / 100.0
    beta = np.linspace(0.0, np.pi, n)
    xc = 0.5 * (1.0 - np.cos(beta))
    yt = 5 * t * (0.2969 * np.sqrt(xc) - 0.1260 * xc - 0.3516 * xc**2 + 0.2843 * xc**3 - 0.1015 * xc**4)
    if m > 0:
        yc = np.where(xc < p, m / p**2 * (2 * p * xc - xc**2), m / (1 - p)**2 * ((1 - 2 * p) + 2 * p * xc - xc**2))
        dyc = np.where(xc < p, 2 * m / p**2 * (p - xc), 2 * m / (1 - p)**2 * (p - xc))
        th = np.arctan(dyc)
    else:
        yc = np.zeros_like(xc)
        th = np.zeros_like(xc)
    xu, yu = xc - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = xc + yt * np.sin(th), yc - yt * np.cos(th)
    return np.concatenate([xu[::-1], xl[1:]]), np.concatenate([yu[::-1], yl[1:]])


def read_dat(path):
    """Reads a Selig-format .dat file (optional name line), returns x, y."""
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.replace(",", " ").split()
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except (ValueError, IndexError):
                continue
    xy = np.array(rows)
    return xy[:, 0], xy[:, 1]


def write_dat(path, x, y, name="XFOILAD airfoil"):
    with open(path, "w") as f:
        f.write(name + "\n")
        for xi, yi in zip(x, y):
            f.write(f"{xi: .8f} {yi: .8f}\n")


def split_surfaces(x, y):
    """Splits counterclockwise coordinates at the LE (min x) into upper and lower
    surfaces, both ordered from LE to TE."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if np.sum(0.5 * (y + np.roll(y, -1)) * (x - np.roll(x, -1))) < 0:   # clockwise -> reverse
        x, y = x[::-1], y[::-1]
    ile = int(np.argmin(x))
    xu, yu = x[:ile + 1][::-1], y[:ile + 1][::-1]
    xl, yl = x[ile:], y[ile:]
    return (xu, yu), (xl, yl)


def normalize(x, y):
    """Translates/rotates/scales so that LE=(0,0) and TE midpoint=(1,0)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ile = int(np.argmin(x))
    xle, yle = x[ile], y[ile]
    xte, yte = 0.5 * (x[0] + x[-1]), 0.5 * (y[0] + y[-1])
    c = np.hypot(xte - xle, yte - yle)
    ang = np.arctan2(yte - yle, xte - xle)
    ca, sa = np.cos(-ang), np.sin(-ang)
    xr, yr = (x - xle), (y - yle)
    return (ca * xr - sa * yr) / c, (sa * xr + ca * yr) / c


def min_thickness(x, y, n=400):
    """Smallest local thickness y_upper(x) - y_lower(x) of a section given by
    its nodes (counterclockwise from the upper TE) and its x/c; negative when
    the lower side crosses above the upper side."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ile = int(np.argmin(x))
    xs = np.linspace(x[ile], min(x[0], x[-1]), n)[1:-1]
    t = np.interp(xs, x[:ile + 1][::-1], y[:ile + 1][::-1]) - np.interp(xs, x[ile:], y[ile:])
    i = int(np.argmin(t))
    return float(t[i]), float(xs[i])

def repanel(x, y, n_side=120):
    """Cubic-spline (arc length) repaneling with a cosine node distribution on
    each side, as XFOIL's PANE.  XFOILAD uses the nodes it is given as panel
    nodes; coordinates coming from elsewhere (other codes, .dat files) may be
    too coarse at the TE: on a NACA 0012 with 220 externally generated nodes, CL at
    4 deg is 6 % low until the section is repaneled."""
    from scipy.interpolate import CubicSpline
    x, y = normalize(x, y)
    s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    sx, sy = CubicSpline(s, x), CubicSpline(s, y)
    s_le = s[int(np.argmin(x))]
    for _ in range(20):                                  # LE = min x of the spline
        s_le -= sx(s_le, 1) / sx(s_le, 2)
    c = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n_side + 1)))
    ss = np.concatenate([s_le * c, s_le + (s[-1] - s_le) * c[1:]])
    return sx(ss), sy(ss)


# --------------------------------------------------------------------------
# base class
# --------------------------------------------------------------------------
class Parametrization:
    """Affine map dv -> panel nodes; subclasses fill x0, y0, Jx, Jy, dv0."""

    names = ()

    @property
    def ndv(self):
        return len(self.dv0)

    @property
    def n_nodes(self):
        return len(self.x0)

    def coordinates(self, dv):
        dv = np.asarray(dv)
        if not np.iscomplexobj(dv):            # complex allowed (complex-step derivatives)
            dv = dv.astype(float)
        return self.x0 + self.Jx @ dv, self.y0 + self.Jy @ dv

    def jacobian(self):
        """(dx/ddv, dy/ddv), each of shape (n_nodes, ndv)."""
        return self.Jx, self.Jy

    # linear operators y_surface(x) = c + row @ dv (to build geometric constraints)
    def upper_op(self, xq):
        raise NotImplementedError

    def lower_op(self, xq):
        raise NotImplementedError

    def thickness_op(self, xq):
        cu, ru = self.upper_op(xq)
        cl, rl = self.lower_op(xq)
        return cu - cl, ru - rl

    def camber_op(self, xq):
        cu, ru = self.upper_op(xq)
        cl, rl = self.lower_op(xq)
        return 0.5 * (cu + cl), 0.5 * (ru + rl)


# --------------------------------------------------------------------------
# B-spline parametrization
# --------------------------------------------------------------------------
def _clamped_knots(n_ctrl, k):
    inner = np.linspace(0.0, 1.0, n_ctrl - k + 1)
    return np.concatenate([np.zeros(k), inner, np.ones(k)])


def _node_spacing(n, r_le=0.08, w_le=0.12, r_te=0.35, w_te=0.08):
    """Normalized arc-length positions of n nodes on one surface (LE->TE),
    clustered at the LE (ratio r_le of the mean panel length) and at the TE."""
    t = np.linspace(0.0, 1.0, 4001)
    g = 1.0 - (1.0 - r_le) * np.exp(-(t / w_le) ** 2) - (1.0 - r_te) * np.exp(-((1.0 - t) / w_te) ** 2)
    s = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(t))])
    s /= s[-1]
    return np.interp(np.linspace(0.0, 1.0, n), np.linspace(0.0, 1.0, len(t)), s)


class BSplineAirfoil(Parametrization):
    """Upper and lower clamped B-splines of degree ``degree`` with ``n_ctrl``
    control points each:

        P0 = (0, 0)          leading edge (shared)
        P1 = (0, +a / -a)    vertical tangent at LE, shared parameter a
        P2..P(n-2)           free y (design variables), fixed x
        P(n-1) = (1, +/- te/2)  trailing edge (fixed thickness ``te``)

    Design vector: [a, yu_2..yu_(n-2), yl_2..yl_(n-2)], or [a, yu_2..yu_(n-2)]
    with ``symmetric=True`` (lower control points yl_i = -yu_i).
    """

    #: max nodes per side: 2*n_side - 1 panel nodes must fit in the DLL (nmax = 364)
    N_SIDE_MAX = 182

    def __init__(self, n_ctrl=8, degree=3, te=0.0025, n_side=None, x_ctrl=None,
                 dv0=None, le_cluster=0.08, te_cluster=0.35, symmetric=False, le_node=None):
        if n_ctrl < degree + 1 or n_ctrl < 4:
            raise ValueError("n_ctrl must be >= max(4, degree+1)")
        self.k = degree
        self.nc = n_ctrl
        self.te = max(float(te), 2e-5)        # exactly sharp TE has a singular derivative
        if n_side is None:                    # ~3 nodes per control point, within the DLL limit
            n_side = min(max(80, 3 * n_ctrl), self.N_SIDE_MAX)
        if n_side > self.N_SIDE_MAX:
            raise ValueError(f"n_side must be <= {self.N_SIDE_MAX} (DLL limit of 364 panel nodes)")
        self.n_side = n_side
        if x_ctrl is None:
            th = np.linspace(0.0, np.pi, n_ctrl - 1)[1:]
            x_ctrl = np.concatenate([[0.0, 0.0], 0.5 * (1.0 - np.cos(th))])
            # keep the first free point away from the LE, but never beyond the next
            # one (with many control points the cosine spacing puts several below 0.02)
            if 0.02 < x_ctrl[3]:
                x_ctrl[2] = max(x_ctrl[2], 0.02)
        self.xc = np.asarray(x_ctrl, float)
        if not (len(self.xc) == n_ctrl and self.xc[0] == 0 and self.xc[1] == 0 and self.xc[-1] == 1):
            raise ValueError("x_ctrl must have n_ctrl values with x0 = x1 = 0 and x(n-1) = 1")
        if np.any(np.diff(self.xc[1:]) <= 0):
            raise ValueError("x_ctrl[1:] must be strictly increasing (otherwise the surface folds back)")
        self.t = _clamped_knots(n_ctrl, degree)
        self.ni = n_ctrl - 3                  # free interior points per side
        self.symmetric = bool(symmetric)
        self.le_node = (not self.symmetric) if le_node is None else bool(le_node)
        self.names = (["a_le"] + [f"yu{i}" for i in range(2, n_ctrl - 1)]
                      + ([] if self.symmetric else [f"yl{i}" for i in range(2, n_ctrl - 1)]))
        self.le_cluster, self.te_cluster = le_cluster, te_cluster
        if dv0 is None:
            dv0 = self._default_dv()
        self.dv0 = np.asarray(dv0, float)
        self._build_nodes()

    # -- control polygon --------------------------------------------------
    def _default_dv(self):
        x, y = naca4("0012", 301)
        return self.fit(x, y, set_dv0=False)

    @classmethod
    def n_dv(cls, n_ctrl, symmetric=False):
        """Number of design variables."""
        return 1 + (1 if symmetric else 2) * (n_ctrl - 3)

    def ctrl_y(self, dv):
        """Control point y for the upper and lower curves."""
        dv = np.asarray(dv, float)
        (c0u, Mu), (c0l, Ml) = self._ctrl_maps()
        return c0u + Mu @ dv, c0l + Ml @ dv

    def _ctrl_maps(self):
        """Affine maps dv -> control y: (c0u, Mu), (c0l, Ml)."""
        nd = self.n_dv(self.nc, self.symmetric)
        Mu = np.zeros((self.nc, nd))
        Ml = np.zeros((self.nc, nd))
        Mu[1, 0], Ml[1, 0] = 1.0, -1.0
        for i in range(self.ni):
            Mu[2 + i, 1 + i] = 1.0
            if self.symmetric:
                Ml[2 + i, 1 + i] = -1.0
            else:
                Ml[2 + i, 1 + self.ni + i] = 1.0
        c0u = np.zeros(self.nc)
        c0l = np.zeros(self.nc)
        c0u[-1], c0l[-1] = 0.5 * self.te, -0.5 * self.te
        return (c0u, Mu), (c0l, Ml)

    def _basis(self, u):
        return BSpline.design_matrix(np.clip(u, 0.0, 1.0), self.t, self.k).toarray()

    def x_of_u(self, u):
        return self._basis(u) @ self.xc

    def u_of_x(self, xq):
        """Inverse of the (monotonic, dv-independent) x(u) map."""
        uf = np.linspace(0.0, 1.0, 20001)
        xf = self.x_of_u(uf)
        return np.interp(np.clip(xq, 0.0, 1.0), xf, uf)

    # -- nodes ----------------------------------------------------------
    def _build_nodes(self):
        """Fixes the node parameters u (from the arc length of the initial
        shape, clustered at LE/TE) and builds the affine node map."""
        (c0u, Mu), (c0l, Ml) = self._ctrl_maps()
        uf = np.linspace(0.0, 1.0, 20001)
        Bf = self._basis(uf)
        xf = Bf @ self.xc
        s_target = _node_spacing(self.n_side, self.le_cluster, 0.12, self.te_cluster, 0.08)
        us = []
        for c0, M in ((c0u, Mu), (c0l, Ml)):
            yf = Bf @ (c0 + M @ self.dv0)
            s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(xf), np.diff(yf)))])
            us.append(np.interp(s_target * s[-1], s, uf))
        if not self.le_node:                  # LE node -> one node per side, half way along the first panel
            us = [np.concatenate([[0.5 * u[1]], u[1:]]) for u in us]
        self.u_upper, self.u_lower = us
        Bu, Bl = self._basis(self.u_upper), self._basis(self.u_lower)
        xu, xl = Bu @ self.xc, Bl @ self.xc
        # upper TE -> LE (reversed), then lower LE -> TE without repeating the LE node
        k = 1 if self.le_node else 0
        self.x0 = np.concatenate([xu[::-1], xl[k:]])
        y0u, Ju = Bu @ c0u, Bu @ Mu
        y0l, Jl = Bl @ c0l, Bl @ Ml
        self.y0 = np.concatenate([y0u[::-1], y0l[k:]])
        self.Jy = np.vstack([Ju[::-1], Jl[k:]])
        self.Jx = np.zeros_like(self.Jy)

    # -- surface operators ---------------------------------------------
    def upper_op(self, xq):
        (c0u, Mu), _ = self._ctrl_maps()
        B = self._basis(self.u_of_x(np.atleast_1d(xq)))
        return B @ c0u, B @ Mu

    def lower_op(self, xq):
        _, (c0l, Ml) = self._ctrl_maps()
        B = self._basis(self.u_of_x(np.atleast_1d(xq)))
        return B @ c0l, B @ Ml

    def le_radius(self, dv):
        """Leading-edge radius of curvature (identical on both sides)."""
        cu, _ = self.ctrl_y(dv)
        bx, by = BSpline(self.t, self.xc, self.k), BSpline(self.t, cu, self.k)
        d1 = np.array([bx(0.0, 1), by(0.0, 1)])
        d2 = np.array([bx(0.0, 2), by(0.0, 2)])
        kappa = abs(d1[0] * d2[1] - d1[1] * d2[0]) / np.linalg.norm(d1) ** 3
        return 1.0 / kappa

    # -- fitting ---------------------------------------------------------
    def fit(self, x, y, set_dv0=True):
        """Least-squares fit of the control points to coordinates (x, y)."""
        x, y = normalize(x, y)
        (xu, yu), (xl, yl) = split_surfaces(x, y)
        self.te = max(float(yu[-1] - yl[-1]), 2e-5) if not hasattr(self, "_te_fixed") else self.te
        (c0u, Mu), (c0l, Ml) = self._ctrl_maps()
        rows, rhs = [], []
        for xs, ys, c0, M in ((xu, yu, c0u, Mu), (xl, yl, c0l, Ml)):
            keep = (xs > 1e-6) & (xs < 1.0 - 1e-6)
            B = self._basis(self.u_of_x(xs[keep]))
            w = 1.0 / np.sqrt(np.maximum(xs[keep], 1e-3))[:, None]   # emphasize the LE
            rows.append(w * (B @ M))
            rhs.append(w[:, 0] * (ys[keep] - B @ c0))
        A, b = np.vstack(rows), np.concatenate(rhs)
        dv, *_ = np.linalg.lstsq(A, b, rcond=None)
        if set_dv0:
            self.dv0 = dv
            self._build_nodes()
        return dv

    @classmethod
    def from_coordinates(cls, x, y, n_ctrl=8, te=None, **kw):
        """Least-squares fit to coordinates (with ``symmetric=True``, the
        symmetric section closest to them, i.e. its thickness distribution)."""
        obj = cls(n_ctrl=n_ctrl, dv0=np.zeros(cls.n_dv(n_ctrl, kw.get("symmetric", False))), **kw)
        if te is not None:
            obj.te = max(float(te), 2e-5)
            obj._te_fixed = True
        obj.fit(x, y, set_dv0=True)
        return obj

    @classmethod
    def from_naca(cls, code="0012", n_ctrl=8, **kw):
        x, y = naca4(code, 301)
        return cls.from_coordinates(x, y, n_ctrl=n_ctrl, **kw)

    @classmethod
    def from_dat(cls, path, n_ctrl=8, **kw):
        x, y = read_dat(path)
        return cls.from_coordinates(x, y, n_ctrl=n_ctrl, **kw)

    def smoothness_op(self, side="both"):
        """Second differences (index-wise) of the control y *perturbation*
        from dv0, P1..P(n-1) on each side: g = D2 @ M @ (dv - dv0).  Linear
        operator (c, R) with g = c + R @ dv."""
        (_, Mu), (_, Ml) = self._ctrl_maps()
        m = self.nc - 1
        D2 = np.zeros((m - 2, m))
        for i in range(m - 2):
            D2[i, i:i + 3] = (1.0, -2.0, 1.0)
        rows = []
        if side in ("both", "upper") or self.symmetric:     # symmetric: the lower side is the mirror
            rows.append(D2 @ Mu[1:])
        if side in ("both", "lower") and not self.symmetric:
            rows.append(D2 @ Ml[1:])
        if not rows:
            raise ValueError("side must be 'both', 'upper' or 'lower'")
        R = np.vstack(rows)
        return -R @ self.dv0, R

    def control_polygon(self, dv):
        cu, cl = self.ctrl_y(dv)
        return (self.xc, cu), (self.xc, cl)


# --------------------------------------------------------------------------
# CST parametrization
# --------------------------------------------------------------------------
class CSTAirfoil(Parametrization):
    """Kulfan CST: y = C(x) * sum_i A_i B_i(x) + x * te/2 (upper), with
    C(x) = sqrt(x) (1 - x), Bernstein polynomials of degree ``order``.
    Design vector: [Au_0..Au_n, Al_0..Al_n], or [Au_0..Au_n] with
    ``symmetric=True`` (Al_i = -Au_i)."""

    def __init__(self, order=6, te=0.0025, n_side=80, dv0=None, symmetric=False, le_node=None):
        self.order = order
        self.symmetric = bool(symmetric)
        self.le_node = (not self.symmetric) if le_node is None else bool(le_node)
        self.te = max(float(te), 2e-5)
        self.n_side = n_side
        beta = np.linspace(0.0, np.pi, n_side)
        self.xs = 0.5 * (1.0 - np.cos(beta))                 # cosine spacing in x
        if not self.le_node:                                 # LE node -> one node per side at x1/2
            self.xs[0] = 0.5 * self.xs[1]
        self.names = [f"Au{i}" for i in range(order + 1)] + (
            [] if self.symmetric else [f"Al{i}" for i in range(order + 1)])
        if dv0 is None:
            x, y = naca4("0012", 301)
            dv0 = self.fit(x, y, set_dv0=False)
        self.dv0 = np.asarray(dv0, float)
        self._build()

    def _shape_rows(self, x):
        from scipy.special import comb
        x = np.atleast_1d(x)
        n = self.order
        C = np.sqrt(x) * (1.0 - x)
        return np.stack([C * comb(n, i) * x**i * (1.0 - x)**(n - i) for i in range(n + 1)], axis=1)

    @classmethod
    def n_dv(cls, order, symmetric=False):
        """Number of design variables."""
        return (1 if symmetric else 2) * (order + 1)

    def _rows(self, xq):
        """Rows of the upper and lower surface maps y = c + row @ dv."""
        S = self._shape_rows(xq)
        if self.symmetric:
            return S, -S
        return np.hstack([S, np.zeros_like(S)]), np.hstack([np.zeros_like(S), S])

    def _build(self):
        Ju, Jl = self._rows(self.xs)
        y0u, y0l = self.xs * 0.5 * self.te, -self.xs * 0.5 * self.te
        k = 1 if self.le_node else 0
        self.x0 = np.concatenate([self.xs[::-1], self.xs[k:]])
        self.y0 = np.concatenate([y0u[::-1], y0l[k:]])
        self.Jy = np.vstack([Ju[::-1], Jl[k:]])
        self.Jx = np.zeros_like(self.Jy)
        assert self.Jy.shape[1] == self.n_dv(self.order, self.symmetric)

    def upper_op(self, xq):
        return np.atleast_1d(xq) * 0.5 * self.te, self._rows(xq)[0]

    def lower_op(self, xq):
        return -np.atleast_1d(xq) * 0.5 * self.te, self._rows(xq)[1]

    def fit(self, x, y, set_dv0=True):
        x, y = normalize(x, y)
        (xu, yu), (xl, yl) = split_surfaces(x, y)
        self.te = max(float(yu[-1] - yl[-1]), 2e-5)
        rows, rhs = [], []
        for xs, ys, op in ((xu, yu, self.upper_op), (xl, yl, self.lower_op)):
            keep = (xs > 1e-6) & (xs < 1 - 1e-6)
            c, R = op(xs[keep])
            rows.append(R)
            rhs.append(ys[keep] - c)
        dv, *_ = np.linalg.lstsq(np.vstack(rows), np.concatenate(rhs), rcond=None)
        if set_dv0:
            self.dv0 = dv
            self._build()
        return dv

    @classmethod
    def from_coordinates(cls, x, y, order=6, **kw):
        obj = cls(order=order, dv0=np.zeros(cls.n_dv(order, kw.get("symmetric", False))), **kw)
        obj.fit(x, y)
        return obj

    @classmethod
    def from_naca(cls, code="0012", order=6, **kw):
        return cls.from_coordinates(*naca4(code, 301), order=order, **kw)

    @classmethod
    def from_dat(cls, path, order=6, **kw):
        return cls.from_coordinates(*read_dat(path), order=order, **kw)


# --------------------------------------------------------------------------
# trailing-edge flap
# --------------------------------------------------------------------------
def _smoothstep(s):
    """Quintic ramp 0 -> 1 on [0, 1] (C2) and its derivative."""
    s = np.clip(s, 0.0, 1.0)
    return s**3 * (10.0 - 15.0 * s + 6.0 * s**2), 30.0 * s**2 * (1.0 - s) ** 2


class FlapDeflection:
    """Plain trailing-edge flap: the nodes aft of the hinge (x_h, y_h) are
    rotated about it by ``delta`` (deg, positive trailing edge down, as
    XFOIL's FLAP).  The rotation angle ramps smoothly from 0 to delta over
    x_h - blend < x < x_h + blend (a faired hinge), so that the deflected
    surface stays smooth and no node overlaps on the compressed side; the
    map is differentiable in delta, the hinge (x_h, y_h) and the nodes.

        theta_i = delta * w(x_i)
        x' = x_h + cos(theta) (x - x_h) + sin(theta) (y - y_h)
        y' = y_h - sin(theta) (x - x_h) + cos(theta) (y - y_h)

    ``pullback`` turns the gradient w.r.t. the deflected nodes (XFOIL
    adjoint) into gradients w.r.t. the undeflected nodes, y_h, delta and x_h."""

    def __init__(self, x, y, delta, x_hinge, y_hinge, blend=0.03):
        if blend <= 0:
            raise ValueError("flap blend must be > 0")
        self.delta, self.xh, self.yh = float(delta), float(x_hinge), float(y_hinge)
        k = np.pi / 180.0
        self.w, wx = _smoothstep((np.asarray(x, float) - self.xh + blend) / (2.0 * blend))
        wx = wx / (2.0 * blend)
        th = k * self.delta * self.w
        self.c, self.s = np.cos(th), np.sin(th)
        dx, dy = x - self.xh, y - self.yh
        self.x = self.xh + self.c * dx + self.s * dy
        self.y = self.yh - self.s * dx + self.c * dy
        thx = k * self.delta * wx
        self.xf_x = self.c + (self.y - self.yh) * thx       # d x'/d x
        self.yf_x = -self.s - (self.x - self.xh) * thx      # d y'/d x
        self.xf_d = k * self.w * (self.y - self.yh)          # d x'/d delta
        self.yf_d = -k * self.w * (self.x - self.xh)         # d y'/d delta

    def pullback(self, gx, gy):
        """(dJ/dx, dJ/dy, dJ/dy_h, dJ/ddelta, dJ/dx_h) from dJ/dx', dJ/dy'.
        The map depends on x and x_h through x - x_h only (rotation and
        ramp), plus x_h itself: dx'/dx_h = 1 - dx'/dx, dy'/dx_h = -dy'/dx."""
        gxs = gx * self.xf_x + gy * self.yf_x
        gys = gx * self.s + gy * self.c
        gyh = float(np.sum(-gx * self.s + gy * (1.0 - self.c)))
        gd = float(np.sum(gx * self.xf_d + gy * self.yf_d))
        gxh = float(np.sum(gx * (1.0 - self.xf_x) - gy * self.yf_x))
        return gxs, gys, gyh, gd, gxh

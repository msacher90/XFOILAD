"""Multi-point airfoil design problem: operating points, objective, hydrodynamic
and geometric constraints, all with adjoint / exact gradients w.r.t. the
shape design variables."""
import os
import warnings
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Optional, Sequence, Union

import numpy as np

from . import geomfuncs
from .geometry import FlapDeflection, min_thickness
from .sizing import FixedChord
from .solver import XFoilAD


class EvaluationError(RuntimeError):
    """XFOIL did not converge for some operating point."""


def ncrit_from_tu(tu):
    """Critical amplification factor from the freestream turbulence level Tu (%),
    Mack's correlation used by XFOIL:  Ncrit = -8.43 - 2.4 ln(Tu/100)."""
    return -8.43 - 2.4 * np.log(np.asarray(tu, float) / 100.0)


def cavitation_number(speed, depth=0.0, p_atm=101325.0, p_v=2340.0, rho=1025.0, g=9.81):
    """Cavitation number sigma = (p_atm + rho g depth - p_v) / (rho speed^2 / 2)
    of a section at ``depth`` (m) below the free surface, at ``speed`` (m/s);
    defaults: sea water at about 20 degC."""
    return (p_atm + rho * g * depth - p_v) / (0.5 * rho * np.asarray(speed, float) ** 2)


def tu_from_ncrit(ncrit):
    """Turbulence level Tu (%) equivalent to Ncrit (inverse of ``ncrit_from_tu``)."""
    return 100.0 * np.exp(-(np.asarray(ncrit, float) + 8.43) / 2.4)


@dataclass
class OperatingPoint:
    """One XFOIL operating point. Give ``cl`` (XFOIL specified-CL mode: alpha
    becomes an output and CD is evaluated at that exact CL) or ``alpha`` (deg).

    Transition: give ``ncrit`` (default 9) or ``tu``, the freestream turbulence
    level in % (converted to Ncrit by Mack's correlation); either may be a
    scalar or a pair (upper, lower).

    Flap: ``flap`` is the deflection (deg, trailing edge down > 0) about the
    hinge of the ``Problem``.  With ``flap_var`` it is a design variable
    (``flap`` is then its initial value): True gives the point its own
    variable, a string names a variable shared by all the points that give
    it (e.g. the samples of a robust design).  At given CL, XFOIL solves for
    alpha and the optimizer chooses the deflection.

    Incidence offset: ``alpha_from`` (index or name of a specified-CL point
    of the same Problem) and ``dalpha`` (deg) make an alpha-mode point at
    alpha = alpha(trim point) + dalpha, e.g. a gust or a wave orbital
    velocity around the trimmed incidence; the trim angle depends on the
    design, and its adjoint gradient is chained into the gradients of
    the offset point.

    ``in_objective=False``: the point does not enter the objective, it only
    carries constraints (e.g. cfmin >= 0 at a specified CL
    and Tu where the drag is not to be minimized).

    ``speed`` (m/s, see ``sizing.KNOT``): with a ``Problem(chord=..., nu=...)``,
    the Reynolds number is Re = speed * chord / nu at each design (``re`` is
    then overwritten), and its dependence on the chord is chained into the
    gradients (exact adjoint derivative w.r.t. Re)."""
    cl: Optional[float] = None
    alpha: Optional[float] = None
    re: float = 1e6
    mach: float = 0.0
    ncrit: Optional[Union[float, Sequence[float]]] = None
    xtr: Sequence[float] = (1.0, 1.0)
    weight: float = 1.0
    name: Optional[str] = None
    niter: int = 150
    tu: Optional[Union[float, Sequence[float]]] = None
    flap: float = 0.0
    flap_var: Optional[Union[bool, str]] = None
    alpha_from: Optional[Union[int, str]] = None
    dalpha: float = 0.0
    in_objective: bool = True
    speed: Optional[float] = None

    def __post_init__(self):
        if self.alpha_from is not None:
            if self.cl is not None:
                raise ValueError("OperatingPoint: alpha_from gives an alpha-mode point, do not give cl")
            self.dalpha = float(self.dalpha)
            if self.alpha is None:
                self.alpha = self.dalpha          # replaced by alpha(trim) + dalpha at each evaluation
            if self.name is None:
                self.name = f"{self.alpha_from}{self.dalpha:+g}deg"
        if (self.cl is None) == (self.alpha is None):
            raise ValueError("OperatingPoint: give exactly one of cl / alpha")
        if self.name is None:
            self.name = f"cl={self.cl}" if self.cl is not None else f"alpha={self.alpha}"
        self.flap = float(self.flap)
        if self.flap_var is True:
            self.flap_var = f"flap@{self.name}"
        elif not self.flap_var:
            self.flap_var = None
        self.xtr = tuple(self.xtr)
        if self.tu is not None:
            if self.ncrit is not None:
                raise ValueError("OperatingPoint: give ncrit or tu, not both")
            tu = np.broadcast_to(np.asarray(self.tu, float), (2,))
            if np.any(tu <= 0):
                raise ValueError("OperatingPoint: tu must be > 0 (in %)")
            nc = ncrit_from_tu(tu)
            if np.any(tu < 0.01) or np.any(nc < 1.0):
                warnings.warn(f"[{self.name}] Tu={self.tu}% gives Ncrit={nc.round(2).tolist()}: "
                              "outside the usual range of Mack's correlation (Tu ~ 0.03 to 2 %)")
            self.ncrit = float(nc[0]) if np.isscalar(self.tu) else (float(nc[0]), float(nc[1]))
        elif self.ncrit is None:
            self.ncrit = 9.0
        if not np.isscalar(self.ncrit):
            self.ncrit = tuple(float(v) for v in self.ncrit)

    @property
    def ncrit_pair(self):
        return (self.ncrit, self.ncrit) if np.isscalar(self.ncrit) else tuple(self.ncrit)

    @property
    def conditions(self):
        """One-line summary of the flow conditions."""
        nt, nb = self.ncrit_pair
        tt, tb = tu_from_ncrit((nt, nb))
        nc = f"Ncrit={nt:.2f}" if nt == nb else f"Ncrit={nt:.2f}/{nb:.2f}"
        tu = f"Tu={tt:.3g}%" if nt == nb else f"Tu={tt:.3g}/{tb:.3g}%"
        s = f"Re={self.re:.3g}  M={self.mach:g}  {nc} ({tu})"
        if self.xtr != (1.0, 1.0):
            s += f"  xtr={self.xtr[0]:g}/{self.xtr[1]:g}"
        if self.alpha_from is not None:
            s += f"  alpha = alpha({self.alpha_from}) {self.dalpha:+g} deg"
        if self.flap_var:
            s += f"  flap free ({self.flap_var})"
        elif self.flap:
            s += f"  flap={self.flap:g} deg"
        if not self.in_objective:
            s += "  (constraints only)"
        return s

    @property
    def mode(self):
        return "cl" if self.cl is not None else "alpha"


HYDRO_FUNCS = ("cl", "cd", "cdf", "cdp", "cm", "alpha", "cd_cl", "ld", "cpmin", "cfmin", "cmh")
#: nodal distributions; "cp:i" / "cf:i" is the value at panel node i (sensitivities)
NODAL_FUNCS = ("cp", "cf")


def _ks_min(cp, mask, rho):
    """Smooth min of cp over mask (= -KS(-cp), a lower bound of the min) and
    its derivative w.r.t. cp."""
    v = -cp[mask]
    m = v.max()
    e = np.exp(rho * (v - m))
    ks = m + np.log(e.sum()) / rho
    w = np.zeros_like(cp)
    w[mask] = e / e.sum()
    return -ks, w


# ----------------------------------------------------------------------
# one operating point: robust primal solve + adjoint passes.  These functions
# run either in-process or in a worker process (see pool.py).
# ----------------------------------------------------------------------
def _setup(solver, p, tol):
    solver.set_flow(re=p.re, mach=p.mach, ncrit=p.ncrit, xtr=p.xtr, verbose=False)
    solver.set_tolerance(tol)


def solve_point(solver, p, x, y, alpha_guess, tol, state=None, warm_first=False):
    """Converges point p from a cold start.  If it fails: warm start from
    ``state`` (``BLState`` of the last convergence of this point, on a
    nearby section), then an alpha continuation.  A point that converges
    cold never depends on the history.  With ``warm_first``, the warm start
    comes first (faster; the result then depends on the history within the
    convergence tolerance).  Leaves the converged state in the library (for
    adjoint).  Returns (result or None, converged by a fallback)."""
    op = dict(cl=p.cl) if p.mode == "cl" else dict(alpha=p.alpha)
    if state is not None and getattr(state, "n", len(x)) != len(x):
        state = None                       # computed on another number of nodes

    def cold():
        _setup(solver, p, tol)
        if p.mode == "cl":
            solver.set_alpha_guess(alpha_guess)
        return solver.analyze(x, y, niter=p.niter, **op)

    def warm():
        _setup(solver, p, tol)
        solver.set_state(state)
        return solver.analyze(x, y, niter=p.niter, warm=True, **op)

    tries = [warm, cold] if warm_first else [cold, warm]
    for k, attempt in enumerate(tries):
        if attempt is warm and state is None:
            continue
        r = attempt()
        if _usable(r):
            return r, k > 0
    r = _continuation(solver, p, x, y, alpha_guess, tol)
    return r, r is not None


def _usable(r):
    """Converged and physical."""
    return r.converged and _physical(r)


def _physical(r):
    """Finite positive drag and, viscous flow, a pressure drag CD - CDf >= 0.
    XFOIL can converge, from a cold start as well as from a fallback, to a
    spurious solution with CD below the skin-friction drag (seen with
    transition at the leading edge and reversed flow); it is rejected."""
    if not np.isfinite(r.cd):
        return False
    if r.cf is None or not np.any(r.cf):                    # inviscid: cf = 0
        return True
    return r.cd > 0.0 and r.cd >= r.cdf


def _continuation(solver, p, x, y, alpha_guess, tol):
    a_target = p.alpha if p.mode == "alpha" else alpha_guess
    starts = (0.0, 0.5 * a_target)
    if p.alpha_from is not None:      # incidence-offset point: start from the (converged) trim angle
        starts = (p.alpha - p.dalpha,) + starts
    for a0 in starts:
        _setup(solver, p, tol)
        r = solver.analyze(x, y, alpha=a0, niter=p.niter)
        if not r.converged:
            continue
        n = max(1, int(np.ceil(abs(a_target - a0) / 0.75)))
        for a in np.linspace(a0, a_target, n + 1)[1:]:
            r = solver.analyze(x, y, alpha=a, niter=p.niter, warm=True)
            if not r.converged:
                break
        else:
            if p.mode == "cl":
                r = solver.analyze(x, y, cl=p.cl, niter=p.niter, warm=True)
            if _usable(r):
                return r
    return None


def hinge_contour(x, xh):
    """Contour of the flap for the hinge moment, built on the undeflected
    nodes x (counterclockwise): upper TE -> upper cut -> hinge -> lower cut
    -> lower TE.  The cuts are the points of the undeflected surface at
    x = xh; each contour point is (1 - f) node[ia] + f node[ib], or the
    hinge itself (is_h).  Constant for a given section x (x does not depend
    on the design variables for BSplineAirfoil and CSTAirfoil) and hinge
    position; df = df/dxh (the two cuts) for the hinge position as a design
    variable."""
    x = np.asarray(x, float)
    ile = int(np.argmin(x))
    if not (x[ile] < xh < min(x[0], x[-1])):
        raise ValueError(f"flap hinge x/c = {xh} must lie between the LE and the TE")
    k = int(np.nonzero(x[:ile + 1] >= xh)[0].max())          # upper: x[k] >= xh > x[k+1]
    m = ile + int(np.nonzero(x[ile:] < xh)[0].max())          # lower: x[m] < xh <= x[m+1]
    ia = list(range(k + 1)) + [k, 0, m] + list(range(m + 1, len(x)))
    ib = list(range(k + 1)) + [k + 1, 0, m + 1] + list(range(m + 1, len(x)))
    f = ([0.0] * (k + 1) + [(x[k] - xh) / (x[k] - x[k + 1]), 0.0, (xh - x[m]) / (x[m + 1] - x[m])]
         + [0.0] * (len(x) - m - 1))
    is_h = np.zeros(len(ia), bool)
    is_h[k + 2] = True
    df = np.zeros(len(ia))
    df[k + 1] = -1.0 / (x[k] - x[k + 1])
    df[k + 3] = 1.0 / (x[m + 1] - x[m])
    # pressure of segment j (points j, j+1) = wa Cp_j + wb Cp_j+1; the inner
    # faces (cut -> hinge -> cut) carry the pressure of their surface cut
    nseg = len(ia) - 1
    wa, wb = np.full(nseg, 0.5), np.full(nseg, 0.5)
    wa[k + 1], wb[k + 1] = 1.0, 0.0
    wa[k + 2], wb[k + 2] = 0.0, 1.0
    return dict(ia=np.array(ia), ib=np.array(ib), f=np.array(f), is_h=is_h, wa=wa, wb=wb, df=df)


def hinge_moment(x, y, cp, xh, yh, ops):
    """Flap hinge moment coefficient (per unit span, / (q c^2)), XFOIL's HMOM:
    sum over the closed flap contour of  Cp_mid (x_mid dx + y_mid dy)  about
    the hinge, counterclockwise positive (x aft, y up): > 0 tends to raise
    the trailing edge.  Returns (value, d/dcp, d/dx, d/dy, d/dyh, d/dxh),
    with x, y the (deflected) nodes the pressures act on; d/dxh: explicit
    dependence on the hinge position (moment arm and flap cuts)."""
    ia, ib, f, h = ops["ia"], ops["ib"], ops["f"], ops["is_h"]
    g = 1.0 - f
    rx = np.where(h, 0.0, g * x[ia] + f * x[ib] - xh)
    ry = np.where(h, 0.0, g * y[ia] + f * y[ib] - yh)
    pp = np.where(h, 0.0, g * cp[ia] + f * cp[ib])
    s = rx**2 + ry**2
    d = 0.5 * np.diff(s)                                   # x_mid dx + y_mid dy, exact
    p = ops["wa"] * pp[:-1] + ops["wb"] * pp[1:]
    val = float(np.sum(p * d))
    dpp = np.zeros_like(pp)                                # d/d(point pressure)
    dpp[:-1] += ops["wa"] * d
    dpp[1:] += ops["wb"] * d
    e = np.zeros_like(s)                                   # d/ds
    e[:-1] -= 0.5 * p
    e[1:] += 0.5 * p
    e[h] = 0.0
    dpp[h] = 0.0
    n = len(x)
    gcp, gx, gy = np.zeros(n), np.zeros(n), np.zeros(n)
    for w, idx in ((g, ia), (f, ib)):
        np.add.at(gcp, idx, w * dpp)
        np.add.at(gx, idx, w * 2.0 * e * rx)
        np.add.at(gy, idx, w * 2.0 * e * ry)
    gyh = float(-np.sum(2.0 * e * ry))
    df = ops.get("df", np.zeros_like(f))
    drx = np.where(h, 0.0, df * (x[ib] - x[ia]) - 1.0)
    dry = np.where(h, 0.0, df * (y[ib] - y[ia]))
    dpx = np.where(h, 0.0, df * (cp[ib] - cp[ia]))
    gxh = float(np.sum(dpp * dpx) + np.sum(2.0 * e * (rx * drx + ry * dry)))
    return val, gcp, gx, gy, gyh, gxh


def cf_mask(r, xrange=(0.0, 1.0), side="both"):
    """Nodes where the skin-friction function cfmin applies: x in xrange, on
    the given side ('both', 'upper', 'lower': split at the LE, min x), the
    two nodes around the stagnation point excluded (Cf -> 0 there)."""
    x, cp = r.x, r.cp
    n = len(x)
    ile = int(np.argmin(x))
    idx = np.arange(n)
    mask = (x >= xrange[0]) & (x <= xrange[1])
    if side == "upper":
        mask &= idx <= ile
    elif side == "lower":
        mask &= idx >= ile
    elif side != "both":
        raise ValueError("cfmin side must be 'both', 'upper' or 'lower'")
    ist = int(np.argmax(cp))
    nb = [j for j in (ist - 1, ist + 1) if 0 <= j < n]
    mask[ist] = False
    mask[max(nb, key=lambda j: cp[j])] = False
    if not mask.any():
        raise ValueError("cfmin: no panel node in the x range / side")
    return mask


def hydro_seed(name, r, opts, hinge=None):
    """Value of an hydrodynamic function, its adjoint seeds (scalars, Cp, Cf)
    and its explicit derivatives (dF/dx, dF/dy, dF/dy_hinge, dF/dx_hinge) or None.
    opts: cpmin_rho, cpmin_xrange, cfmin_rho, cfmin_xrange, cfmin_side."""
    scal = {"cl": {"cl": 1.0}, "cd": {"cd": 1.0}, "cdf": {"cdf": 1.0}, "cm": {"cm": 1.0},
            "alpha": {"alpha": 1.0}, "cdp": {"cd": 1.0, "cdf": -1.0}}
    if name in scal:
        v = r.cd - r.cdf if name == "cdp" else getattr(r, name)
        return v, scal[name], None, None, None
    if name == "cd_cl":
        return r.cd / r.cl, {"cd": 1.0 / r.cl, "cl": -r.cd / r.cl**2}, None, None, None
    if name == "ld":
        return r.cl / r.cd, {"cl": 1.0 / r.cd, "cd": -r.cl / r.cd**2}, None, None, None
    if name == "cpmin":
        xr = opts["cpmin_xrange"]
        mask = (r.x >= xr[0]) & (r.x <= xr[1])
        v, w = _ks_min(r.cp, mask, opts["cpmin_rho"])
        return v, {}, w, None, None
    if name == "cfmin":
        v, w = _ks_min(r.cf, cf_mask(r, opts["cfmin_xrange"], opts["cfmin_side"]), opts["cfmin_rho"])
        return v, {}, None, w, None
    if name == "cmh":
        if hinge is None:
            raise ValueError("cmh needs the flap hinge of the Problem")
        v, gcp, gx, gy, gyh, gxh = hinge_moment(r.x, r.y, r.cp, *hinge)
        return v, {}, gcp, None, (gx, gy, gyh, gxh)
    base, _, i = name.partition(":")
    if base in NODAL_FUNCS and i.isdigit():
        e = np.zeros(len(r.x))
        e[int(i)] = 1.0
        v = float(getattr(r, base)[int(i)])
        return (v, {}, e, None, None) if base == "cp" else (v, {}, None, e, None)
    raise ValueError(f"unknown hydrodynamic function {name!r}")


def run_point(solver, task):
    """Solves one operating point and runs one adjoint pass per requested
    function.  task: dict(point, x, y, names, gradients, alpha_guess, tol,
    kiter, cf_kiter, opts, hinge, converged).  ``converged`` = (result, BL
    state) of this point already converged at this design (evaluation
    without gradients): the solve is skipped and the adjoint passes restart
    from that state (same gradients, bit for bit).  Returns dict(result,
    values, grads, continuation, error, state); result is None when XFOIL
    failed; grads[name] = (dF/dx, dF/dy, dF/dop, dF/dy_hinge, dF/dRe, dF/dx_hinge)."""
    p, x, y = task["point"], task["x"], task["y"]
    if task.get("converged") is not None and task["gradients"]:
        r, state = task["converged"]
        _setup(solver, p, task["tol"])
        solver.set_state(state)
        cont = False
    else:
        r, cont = solve_point(solver, p, x, y, task["alpha_guess"], task["tol"], task.get("state"),
                              task.get("warm_first", False))
        state = r.state if r is not None else None
    out = dict(result=r, values={}, grads={}, continuation=cont, error=None, state=state)
    if r is None:
        out["error"] = f"XFOIL failed at operating point {p.name}"
        return out
    op = dict(cl=p.cl) if p.mode == "cl" else dict(alpha=p.alpha)
    for name in task["names"]:
        v, seed, cps, cfs, explicit = hydro_seed(name, r, task["opts"], task.get("hinge"))
        out["values"][name] = v
        if task["gradients"]:
            # the skin friction needs more replayed Newton iterations than the
            # integrated outputs and Cp to reach the converged adjoint
            k = max(task["kiter"], task["cf_kiter"]) if cfs is not None else task["kiter"]
            gx, gy, gop, gfl = solver.adjoint(x, y, seed=seed, cp_seed=cps, cf_seed=cfs, kiter=k, flow=True, **op)
            if not (np.all(np.isfinite(gx)) and np.all(np.isfinite(gy))):
                out["error"] = f"non-finite adjoint at {p.name}"
                return out
            gyh = gxh = 0.0
            if explicit is not None:          # direct dependence on the nodes and the hinge
                gx, gy, gyh, gxh = gx + explicit[0], gy + explicit[1], explicit[2], explicit[3]
            out["grads"][name] = (gx, gy, gop, gyh, gfl["re"], gxh)
    return out


class HydroEvaluator:
    """Robust primal solves + adjoint gradients for a list of operating points.

    ``n_workers``: number of worker processes, each with its own copy of the
    XFOIL library.  With gradients, every (point, function) pair is an
    independent task (primal solve + one adjoint pass), so the tasks of one
    design run in parallel; results are identical to the serial ones.
    None = environment variable XFOILAD_WORKERS if set, else half the logical
    CPUs; 1 = serial, in-process.

    ``warm_start``: each point first restarts from its BL state converged at
    the previous design (cold start and continuation as fallbacks), instead
    of a cold start first; see ``Problem``.

    ``cpmin``: smooth min of Cp over cpmin_xrange (KS, cpmin_rho).  ``cfmin``:
    smooth min of the skin friction over cfmin_xrange and cfmin_side (KS,
    cfmin_rho, in 1/Cf: the value lies below the true min by at most
    ln(n)/rho); its adjoint replays cf_kiter Newton iterations."""

    def __init__(self, points, solver=None, kiter=4, tol=1e-5, cpmin_rho=80.0,
                 cpmin_xrange=(0.0, 0.98), verbose=False, n_workers=None, cfmin_rho=5e4,
                 cfmin_xrange=(0.0, 1.0), cfmin_side="both", cf_kiter=6, warm_start=False):
        self.points = list(points)
        self.warm_start = bool(warm_start)
        self._solver = solver
        self.kiter = kiter
        self.tol = tol
        self.cpmin_rho = cpmin_rho
        self.cpmin_xrange = cpmin_xrange
        self.cfmin_rho = cfmin_rho
        self.cfmin_xrange = cfmin_xrange
        self.cfmin_side = cfmin_side
        self.cf_kiter = cf_kiter
        self.verbose = verbose
        if n_workers is None:
            n_workers = os.environ.get("XFOILAD_WORKERS") or max(1, (os.cpu_count() or 2) // 2)
        self.n_workers = int(n_workers)
        self._pool = None
        self._alpha_guess = {i: (p.alpha if p.alpha is not None else 0.0) for i, p in enumerate(self.points)}
        self._state = {}          # BL state of the last convergence of each point (warm-start fallback)

    @property
    def solver(self):
        """In-process solver (serial mode), created on first use."""
        if self._solver is None:
            self._solver = XFoilAD()
        return self._solver

    def close(self):
        """Stops the worker processes (also done automatically at exit)."""
        if self._pool is not None:
            self._pool.close()
            self._pool = None

    def _run(self, tasks):
        if self.n_workers <= 1 or len(tasks) <= 1:
            return [run_point(self.solver, t) for t in tasks]
        if self._pool is None:
            from .pool import WorkerPool
            self._pool = WorkerPool(self.n_workers)
        return self._pool.map(tasks)

    def evaluate(self, x, y, requests, gradients=True, alpha_guess=None, coords=None, hinge=None,
                 converged=None):
        """requests: {point index: set of function names}; alpha_guess:
        optional {point index: alpha} overriding the last converged angles;
        coords: optional {point index: (x, y)} nodes of a point that differ
        from (x, y), e.g. with a deflected flap; hinge: (x_h, y_h,
        hinge_contour) for the hinge moment 'cmh'; converged: optional
        {point index: (result, BL state)} of the points already converged on
        these nodes (gradients only: no new solve, see ``run_point``).
        Returns (results, values, grads) with values[(ip, name)] and
        grads[(ip, name)] = (dF/dx_nodes, dF/dy_nodes, dF/dop, dF/dy_hinge);
        the BL states of the converged points are left in ``last_states``."""
        x, y = np.asarray(x, float), np.asarray(y, float)
        coords = coords or {}
        opts = dict(cpmin_rho=self.cpmin_rho, cpmin_xrange=tuple(self.cpmin_xrange), cfmin_rho=self.cfmin_rho,
                    cfmin_xrange=tuple(self.cfmin_xrange), cfmin_side=self.cfmin_side)
        common = dict(gradients=gradients, tol=self.tol, kiter=self.kiter, cf_kiter=self.cf_kiter,
                      hinge=hinge, opts=opts, warm_first=self.warm_start)
        # in parallel, one task per adjoint pass (each worker re-solves the cheap
        # primal); in serial, one task per point (one primal for all its passes)
        split = gradients and self.n_workers > 1
        tasks, owners = [], []
        for ip in sorted(requests):
            names = sorted(requests[ip])
            xp, yp = coords.get(ip, (x, y))
            for group in ([n] for n in names) if split and names else [names]:
                tasks.append(dict(common, point=self.points[ip], names=group,
                                  x=np.asarray(xp, float), y=np.asarray(yp, float),
                                  alpha_guess=(alpha_guess or self._alpha_guess)[ip],
                                  state=self._state.get(ip), converged=(converged or {}).get(ip)))
                owners.append(ip)
        outs = self._run(tasks)
        results, values, grads, cont = {}, {}, {}, set()
        for ip, out in zip(owners, outs):
            if out["error"]:
                raise EvaluationError(out["error"])
            results.setdefault(ip, out["result"])
            if out["continuation"]:
                cont.add(ip)
            for name, v in out["values"].items():
                values[(ip, name)] = v
            for name, g in out["grads"].items():
                grads[(ip, name)] = g
        for ip, r in results.items():
            self._alpha_guess[ip] = r.alpha
        self.last_states = {}
        for ip, out in zip(owners, outs):
            if out["state"] is not None:
                self._state[ip] = out["state"]
                self.last_states[ip] = out["state"]
            if self.verbose and ip in cont:
                print(f"    [{self.points[ip].name}] converged by warm start or alpha continuation")
        return results, values, grads


@dataclass
class Constraint:
    kind: str
    lower: Optional[float] = None
    upper: Optional[float] = None
    point: Optional[int] = None
    options: dict = field(default_factory=dict)
    scale: float = 1.0
    name: str = ""

    @property
    def is_hydro(self):
        return self.kind in HYDRO_FUNCS


class Problem:
    """Shape optimization problem.

    objective: name of an hydrodynamic function (default 'cd'), or 'drag'
    = chord * CD (drag per unit span and dynamic pressure, m; CD without a
    chord model), combined over the operating points as
    sum_i weight_i * f_i / |f_i(initial)| / sum_i weight_i (normalize='point'),
    or sum_i weight_i * f_i / |sum_i weight_i f_i(initial)| (normalize='global':
    a weighted sum, e.g. a mean over samples of an uncertain parameter); the
    points with ``in_objective=False`` only carry constraints.

    chord: dimensional chord (m), a number or a model of ``sizing``
    (``ChordFromInertia``, ``ChordFromModulus``: the chord depends on the
    shape); nu: kinematic viscosity (m^2/s) for the points given by their
    speed (Re = speed * chord / nu).

    warm_start: each operating point first restarts from its BL state
    converged at the previous design (``BLState``), then from scratch and by
    continuation if that fails.  Faster (a few Newton iterations instead of
    about ten), but the solution then depends on the history of the
    optimization within the convergence tolerance ``tol``; the gradients stay
    those of the converged solution.  Default False: a point that converges
    from scratch never depends on the history.

    cfmin_rho, cfmin_xrange, cfmin_side, cf_kiter: settings of the
    skin-friction function cfmin (see ``HydroEvaluator``).

    Design vector: the shape variables of ``param``, then one flap deflection
    (deg) per flap variable named by the points (``OperatingPoint.flap_var``).
    The flap hinge is at x/c = ``flap_hinge[0]``, at the fraction
    ``flap_hinge[1]`` of the local thickness above the lower surface (XFOIL's
    y/t; it moves with the shape).  ``flap_bounds`` = (lower, upper) in deg,
    or a dict {variable name: (lower, upper)}; ``flap_blend`` is the half
    width of the faired hinge (see ``geometry.FlapDeflection``);
    ``flap_scale`` is the optimizer scaling of the deflections (deg).
    ``hinge_var=True`` makes the hinge position x/c a design variable (last
    of the design vector, named 'hinge_x', initial value ``flap_hinge[0]``,
    bounds ``hinge_bounds``, optimizer scaling ``hinge_scale``), shared by all
    the points with a flap; its gradient is exact (rotation, hinge height and
    hinge-moment contour), the slope of the surfaces at the hinge being taken
    by a central difference of step 1e-6.
    Geometric constraints apply to the undeflected section."""

    CACHE_SIZE = 16

    def __init__(self, param, points, objective="cd", dv_bounds=0.04, kiter=4, tol=1e-5,
                 solver=None, cpmin_rho=80.0, cpmin_xrange=(0.0, 0.98), verbose=True, n_workers=None,
                 flap_hinge=(0.75, 0.5), flap_bounds=(-10.0, 15.0), flap_blend=0.03, flap_scale=1.0,
                 hinge_var=False, hinge_bounds=(0.55, 0.85), hinge_scale=0.01,
                 cfmin_rho=5e4, cfmin_xrange=(0.0, 1.0), cfmin_side="both", cf_kiter=6,
                 chord=None, nu=None, normalize="point", warm_start=False):
        self.param = param
        self.points = [p if isinstance(p, OperatingPoint) else OperatingPoint(**p) for p in points]
        if objective not in HYDRO_FUNCS + ("drag",):
            raise ValueError(f"objective must be one of {HYDRO_FUNCS + ('drag',)}")
        if normalize not in ("point", "global"):
            raise ValueError("normalize must be 'point' or 'global'")
        self.objective, self.normalize = objective, normalize
        self.fobj = "cd" if objective == "drag" else objective      # function solved at each point
        self.chord_model = FixedChord(chord) if np.isscalar(chord) else chord
        self.nu = None if nu is None else float(nu)
        if any(p.speed is not None for p in self.points):
            if self.chord_model is None or self.nu is None:
                raise ValueError("points given by their speed need Problem(chord=..., nu=...)")
        if not self.objective_points:
            raise ValueError("no operating point in the objective (all have in_objective=False)")
        self.hydro = HydroEvaluator(self.points, solver, kiter, tol, cpmin_rho, cpmin_xrange, verbose, n_workers,
                                    cfmin_rho=cfmin_rho, cfmin_xrange=cfmin_xrange, cfmin_side=cfmin_side,
                                    cf_kiter=cf_kiter, warm_start=warm_start)
        self.constraints = []
        self.verbose = verbose
        self.n_shape = param.ndv
        # flap variables, in order of first appearance; initial value from that point
        self.flap_vars, flap0 = [], []
        for p in self.points:
            if p.flap_var and p.flap_var not in self.flap_vars:
                self.flap_vars.append(p.flap_var)
                flap0.append(p.flap)
        self._flap0 = np.array(flap0, float)
        # incidence-offset points: {index: index of their specified-CL trim point}
        names = [p.name for p in self.points]
        self.trim_of = {}
        for ip, p in enumerate(self.points):
            if p.alpha_from is None:
                continue
            it = p.alpha_from if isinstance(p.alpha_from, (int, np.integer)) else (
                names.index(p.alpha_from) if p.alpha_from in names else None)
            if it is None or not 0 <= it < len(self.points):
                raise ValueError(f"[{p.name}] alpha_from={p.alpha_from!r}: no such operating point")
            if self.points[it].mode != "cl" or self.points[it].alpha_from is not None:
                raise ValueError(f"[{p.name}] alpha_from must name a specified-CL point")
            self.trim_of[ip] = int(it)
        self.flap_index = {ip: self.n_shape + self.flap_vars.index(p.flap_var)
                           for ip, p in enumerate(self.points) if p.flap_var}
        self.flap_hinge = (float(flap_hinge[0]), float(flap_hinge[1]))
        self.flap_blend = float(flap_blend)
        self.flap_scale = float(flap_scale)
        self.hinge_var = bool(hinge_var)
        self.hinge_bounds = (float(hinge_bounds[0]), float(hinge_bounds[1]))
        self.hinge_scale = float(hinge_scale)
        if self.hinge_var:
            if not self.has_flap:
                raise ValueError("hinge_var needs at least one point with a flap")
            if not self.hinge_bounds[0] <= self.flap_hinge[0] <= self.hinge_bounds[1]:
                raise ValueError(f"flap_hinge x/c = {self.flap_hinge[0]} outside hinge_bounds {self.hinge_bounds}")
        self.hinge_index = self.n_shape + len(self.flap_vars) if self.hinge_var else None
        self._hinge_cache = {}
        self._hinge_now = None
        self._hinge_geometry(self.flap_hinge[0])         # checks that the hinge lies on the section
        dv0 = param.dv0
        if np.isscalar(dv_bounds):
            lb, ub = dv0 - dv_bounds, dv0 + dv_bounds
        else:
            lb, ub = (np.array(b, float)[:self.n_shape] for b in dv_bounds)
        if hasattr(param, "le_radius"):            # keep a positive LE radius
            lb[0] = max(lb[0], 0.25 * dv0[0])
        fb = [flap_bounds.get(v, (-10.0, 15.0)) if isinstance(flap_bounds, dict) else flap_bounds
              for v in self.flap_vars]
        self.lb = np.concatenate([lb, [b[0] for b in fb], [self.hinge_bounds[0]] if self.hinge_var else []])
        self.ub = np.concatenate([ub, [b[1] for b in fb], [self.hinge_bounds[1]] if self.hinge_var else []])
        self._cache = OrderedDict()
        self.ref = None
        self._chord_now = (1.0, np.zeros(self.ndv))

    # ------------------------------------------------------ design vector
    @property
    def dv0(self):
        return np.concatenate([self.param.dv0, self._flap0, [self.flap_hinge[0]] if self.hinge_var else []])

    @property
    def ndv(self):
        return self.n_shape + len(self.flap_vars) + int(self.hinge_var)

    @property
    def names(self):
        return (list(getattr(self.param, "names", [f"dv{k}" for k in range(self.n_shape)])) + self.flap_vars
                + (["hinge_x"] if self.hinge_var else []))

    @property
    def objective_points(self):
        """Indices of the operating points in the objective."""
        return [ip for ip, p in enumerate(self.points) if p.in_objective]

    @property
    def has_flap(self):
        return any(p.flap_var or p.flap for p in self.points)

    def dv_scales(self, dv_scale):
        """Optimizer scaling of each design variable: dv_scale for the shape,
        flap_scale (deg) for the deflections, hinge_scale for the hinge position."""
        return np.concatenate([np.broadcast_to(np.asarray(dv_scale, float), (self.n_shape,)),
                               np.full(len(self.flap_vars), self.flap_scale),
                               [self.hinge_scale] if self.hinge_var else []])

    def shape_dv(self, dv):
        return np.asarray(dv, float)[:self.n_shape]

    def chord(self, dv):
        """Dimensional chord (1 without chord model) and its gradient w.r.t. dv."""
        g = np.zeros(self.ndv)
        if self.chord_model is None:
            return 1.0, g
        c, dc = self.chord_model(self.param, self.shape_dv(dv))
        g[:self.n_shape] = dc
        return float(c), g

    def flap_angles(self, dv):
        """Flap deflection (deg) of each operating point."""
        dv = np.asarray(dv, float)
        return [float(dv[self.flap_index[ip]]) if ip in self.flap_index else p.flap
                for ip, p in enumerate(self.points)]

    def hinge_x(self, dv):
        """Position x/c of the flap hinge (a design variable with hinge_var)."""
        return float(np.asarray(dv, float)[self.hinge_index]) if self.hinge_var else self.flap_hinge[0]

    def _hinge_geometry(self, xh):
        """At hinge position xh: hinge height y_h = c + row @ d (linear in the
        shape d), its derivative w.r.t. xh, dc + drow @ d (slope of the
        surfaces, central difference of step 1e-6), and the flap contour of
        the hinge moment.  Cached per position."""
        hit = self._hinge_cache.get(xh)
        if hit is None:
            fh = self.flap_hinge[1]

            def height(xq):
                cu, ru = self.param.upper_op(np.array([xq]))
                cl, rl = self.param.lower_op(np.array([xq]))
                return float((1 - fh) * cl[0] + fh * cu[0]), (1 - fh) * rl[0] + fh * ru[0]

            c, row = height(xh)
            h = 1e-6
            (cp_, rp), (cm_, rm) = height(xh + h), height(xh - h)
            ops = hinge_contour(self.param.coordinates(self.param.dv0)[0], xh)
            hit = (c, row, (cp_ - cm_) / (2 * h), (rp - rm) / (2 * h), ops)
            if len(self._hinge_cache) > 64:
                self._hinge_cache.clear()
            self._hinge_cache[xh] = hit
        return hit

    def hinge_y(self, dv):
        """Height of the flap hinge (it moves with the shape)."""
        c, row = self._hinge_geometry(self.hinge_x(dv))[:2]
        return float(c + row @ self.shape_dv(dv))

    def hinge_moments(self, ev):
        """Hinge moment coefficient (XFOIL's HMOM) of each point of an evaluation."""
        xh = self.hinge_x(ev["dv"])
        yh, ops = self.hinge_y(ev["dv"]), self._hinge_geometry(xh)[4]
        return [hinge_moment(r.x, r.y, r.cp, xh, yh, ops)[0] for _, r in sorted(ev["results"].items())]

    def _deflect(self, x, y, d, delta, xh):
        c, row = self._hinge_geometry(xh)[:2]
        fl = FlapDeflection(x, y, delta, xh, c + row @ d, self.flap_blend)
        ile = int(np.argmin(x))           # the deflected sides must stay monotonic in x
        # (lower side from the node after the LE: without LE node, the first two nodes have the same x)
        if np.any(np.diff(fl.x[:ile + 1]) >= 0) or np.any(np.diff(fl.x[ile + 1:]) <= 0):
            raise EvaluationError(f"flap deflection {delta:.2f} deg folds the surface "
                                  f"(increase flap_blend or reduce the flap bounds)")
        return fl

    def coordinates(self, dv, point=None):
        """Nodes of the section (undeflected), or of the section seen by
        operating point ``point`` (flap deflected)."""
        d = self.shape_dv(dv)
        x, y = self.param.coordinates(d)
        if point is None or not (self.points[point].flap_var or self.points[point].flap):
            return x, y
        fl = self._deflect(x, y, d, self.flap_angles(dv)[point], self.hinge_x(dv))
        return fl.x, fl.y

    # ------------------------------------------------------- constraints
    def add_constraint(self, kind, lower=None, upper=None, equals=None, point=None, scale=None,
                       name=None, **options):
        """Adds a constraint ``lower <= f <= upper`` (or ``f == equals``).

        Hydrodynamic kinds (need ``point``, an index, a name or 'all'): cl, cd,
        cdf, cdp, cm, alpha, cd_cl, ld, cpmin, cfmin (smooth min of the skin
        friction: ``lower=0`` forbids separated flow, see ``cf_mask``), cmh
        (flap hinge moment, see ``hinge_moment``).
        Geometric kinds: thickness(x=...), camber(x=...), max_thickness,
        area, inertia(axis='xx'|'yy'|'xy'|'p'), section_modulus,
        box(x_range=(x1,x2), y_range=(ybot,ytop)), le_radius, te_angle,
        smoothness(side='both'|'upper'|'lower') with lower=-s, upper=s."""
        if equals is not None:
            lower = upper = equals
        if kind == "box":
            lower, upper = 0.0, None
        if lower is None and upper is None:
            raise ValueError("constraint without bound")
        if kind in HYDRO_FUNCS:
            names = [p.name for p in self.points]
            if point in (None, "all"):
                pts = range(len(self.points))
            elif isinstance(point, str):
                if point not in names:
                    raise ValueError(f"no operating point named {point!r}")
                pts = [names.index(point)]
            else:
                pts = [int(point)]
            for ip in pts:
                ref = lower if lower is not None else upper
                sc = scale or max(abs(ref), 1e-3)
                self.constraints.append(Constraint(kind, lower, upper, ip, options, sc,
                                                   name or f"{kind}@{self.points[ip].name}"))
        elif kind in geomfuncs.REGISTRY:
            if kind == "box":
                lower, upper = 0.0, None
            ref = lower if lower is not None else upper
            sc = scale or max(abs(ref), 1e-4) if kind != "box" else (scale or 0.01)
            self.constraints.append(Constraint(kind, lower, upper, None, options, sc, name or kind))
        else:
            raise ValueError(f"unknown constraint kind {kind!r}")
        self._cache.clear()
        return self

    # ------------------------------------------------------- evaluation
    def _requests(self):
        req = {ip: ({self.fobj} if p.in_objective else set()) for ip, p in enumerate(self.points)}
        for c in self.constraints:
            if c.is_hydro:
                req[c.point].add(c.kind)
        for it in set(self.trim_of.values()):     # d alpha_trim / d dv, chained into the offset points
            req[it].add("alpha")
        return req

    def _solve(self, dv, req, gradients, guess=None, converged=None):
        """XFOIL solves (and adjoint passes) of the requests {point: names}
        for design dv: trim points first, then the incidence-offset points at
        alpha(trim) + dalpha; converged: {point: (result, BL state)} of an
        evaluation of dv without gradients (no new solve).  Returns (x, y,
        results, values, grads, defl, states)."""
        d = self.shape_dv(dv)
        x, y = self.param.coordinates(d)
        tmin, xmin = min_thickness(x, y)
        if tmin < 0.0:                     # lower side above the upper side: XFOIL results meaningless
            raise EvaluationError(f"the section crosses itself (thickness {tmin:.2e} at x/c = {xmin:.3f})")
        self._chord_now = self.chord(dv)
        for p in self.points:              # Re = V c / nu at this design
            if p.speed is not None:
                p.re = p.speed * self._chord_now[0] / self.nu
        flaps = self.flap_angles(dv)
        xh = self.hinge_x(dv)
        c, row, dc, drow, ops = self._hinge_geometry(xh)
        self._hinge_now = (row, float(dc + drow @ d))    # d y_h / d shape, d y_h / d x_h
        defl = {ip: self._deflect(x, y, d, flaps[ip], xh) for ip, p in enumerate(self.points)
                if p.flap_var or p.flap}
        coords = {ip: (fl.x, fl.y) for ip, fl in defl.items()}
        hinge = (xh, float(c + row @ d), ops)
        stages = [{ip: r for ip, r in req.items() if ip not in self.trim_of},
                  {ip: r for ip, r in req.items() if ip in self.trim_of}]
        results, vals, grads, states = {}, {}, {}, {}
        for st in stages:
            if not st:
                continue
            for ip in st:
                if ip in self.trim_of:
                    p = self.points[ip]
                    p.alpha = results[self.trim_of[ip]].alpha + p.dalpha
            r_, v_, g_ = self.hydro.evaluate(x, y, st, gradients, alpha_guess=guess,
                                             coords=coords, hinge=hinge, converged=converged)
            results.update(r_)
            vals.update(v_)
            grads.update(g_)
            states.update(self.hydro.last_states)
        return x, y, results, vals, grads, defl, states

    def _chain(self, grads, defl, k):
        """Gradient w.r.t. the design vector of function k = (point, name),
        from its gradient w.r.t. the nodes, the operating parameter and the
        hinge height and position (through the flap rotation and the trim
        angle), and the Reynolds number when it follows the chord (Re = V c / nu)."""
        ip = k[0]
        gx, gy, gop, gyh, gre, gxh = grads[k]    # gyh, gxh: direct dependence on the hinge (cmh)
        g = np.zeros(self.ndv)
        if ip in defl:            # through the flap rotation (nodes, hinge, deflection)
            gx, gy, gyh_rot, gd, gxh_rot = defl[ip].pullback(gx, gy)
            gyh += gyh_rot
            gxh += gxh_rot
            if ip in self.flap_index:
                g[self.flap_index[ip]] = gd
        Jx, Jy = self.param.jacobian()
        row, dyh_dxh = self._hinge_now
        g[:self.n_shape] = Jx.T @ gx + Jy.T @ gy + gyh * row
        if self.hinge_var:        # x_h moves the rotation, the contour and the hinge height
            g[self.hinge_index] = gxh + gyh * dyh_dxh
        p = self.points[ip]
        if p.speed is not None:   # Re = V c(dv) / nu
            c, dc = self._chord_now
            g += gre * (p.re / c) * dc
        if ip in self.trim_of:    # alpha = alpha_trim(dv) + dalpha
            g += gop * self._chain(grads, defl, (self.trim_of[ip], "alpha"))
        return g

    def evaluate(self, dv, gradients=True):
        """Returns dict with objective, constraint values (raw) and gradients
        w.r.t. dv; results are cached per design vector."""
        dv = np.asarray(dv, float)
        key = dv.tobytes()
        hit = self._cache.get(key)
        if isinstance(hit, EvaluationError):      # known failure: do not rerun the whole XFOIL cascade
            raise hit
        if hit is not None and (hit["grad"] is not None or not gradients):
            self._cache.move_to_end(key)
            return hit
        d = self.shape_dv(dv)
        flaps = self.flap_angles(dv)
        # design already converged without gradients: the adjoint passes restart
        # from its converged BL states (no new solve, same solution)
        guess, conv = None, None
        if hit is not None:
            guess = {ip: r.alpha for ip, r in hit["results"].items()}
            conv = {ip: (hit["results"][ip], st) for ip, st in hit.get("states", {}).items()}
        try:
            x, y, results, vals, grads, defl, states = self._solve(dv, self._requests(), gradients, guess, conv)
        except EvaluationError as e:
            self._store(key, e)
            raise

        def chain(k):
            return self._chain(grads, defl, k)

        obj = self.objective_points
        chord, dchord = self._chord_now
        # objective function of each point: f_k (chord * CD for 'drag') and its gradient
        fk = np.array([vals[(ip, self.fobj)] * (chord if self.objective == "drag" else 1.0) for ip in obj])
        fk_grad = None
        if gradients:
            fk_grad = np.vstack([chord * chain((ip, "cd")) + vals[(ip, "cd")] * dchord if self.objective == "drag"
                                 else chain((ip, self.fobj)) for ip in obj])
        w = np.array([self.points[ip].weight for ip in obj])
        if self.ref is None:          # normalization of the objective by the initial values
            self.ref = np.abs(fk) if self.normalize == "point" else abs(float(w @ fk))
        if self.normalize == "point":
            f = float(np.sum(w * fk / self.ref) / w.sum())
            df = (w / self.ref / w.sum()) @ fk_grad if gradients else None
        else:
            f = float(w @ fk) / self.ref
            df = (w @ fk_grad) / self.ref if gradients else None

        nflap = self.ndv - self.n_shape           # deflections and hinge position
        cvals, cgrads = [], []
        for c in self.constraints:
            if c.is_hydro:
                v = np.array([vals[(c.point, c.kind)]])
                J = chain((c.point, c.kind))[None, :] if gradients else None
            else:
                v, J = geomfuncs.REGISTRY[c.kind](self.param, d, **c.options)
                if nflap:             # the flap does not change the section
                    J = np.hstack([np.atleast_2d(J), np.zeros((np.atleast_2d(J).shape[0], nflap))])
            cvals.append(np.atleast_1d(v))
            cgrads.append(J)
        out = dict(dv=dv, x=x, y=y, f=f, grad=df, results=results, values=vals, flaps=flaps,
                   cvals=cvals, cgrads=cgrads if gradients else None, chord=chord, fk=fk, fk_grad=fk_grad,
                   re=[p.re for p in self.points], states=states)
        if not gradients:
            out["grad"] = None
        self._store(key, out)
        return out

    def _store(self, key, entry):
        """LRU cache of the last designs (results or EvaluationError): SLSQP
        asks for the objective, each constraint and their gradients separately."""
        self._cache[key] = entry
        self._cache.move_to_end(key)
        while len(self._cache) > self.CACHE_SIZE:
            self._cache.popitem(last=False)

    # ------------------------------------------------------ sensitivities
    def sensitivities(self, dv=None, functions=("cl", "cd", "cp"), points=None, method="adjoint"):
        """Sensitivities of hydrodynamic functions of design ``dv`` (default:
        the initial design), at the operating points ``points`` (indices or
        names, default all).

        ``functions``: scalar functions (names of HYDRO_FUNCS) and/or the
        nodal distributions 'cp' and 'cf'.

        Returns one dict per point: name, mode ('alpha' or 'cl'), x, y (nodes
        seen by the point, flap deflected), values {f: value or nodal array},
        dx, dy {f: dF/dx_j, dF/dy_j w.r.t. those nodes, at the fixed operating
        parameter of the point; for 'cp' and 'cf' the (n, n) matrices
        [i, j] = dF_i/dx_j}, dop {f: dF/dalpha (deg) in alpha mode, dF/dCL in
        CL mode}, ddv {f: dF/d(design vector), full chain: shape, flap, trim
        angle of the incidence-offset points, Re through the chord}.

        method='adjoint': one adjoint pass per scalar function and, for 'cp'
        and 'cf', per panel node.  method='tangent': the nodal distributions
        are differentiated by the tangent library (``XFoilTangent``), one pass
        per design variable for all of them (ddv only: their dx, dy, dop are
        None); the scalar functions still use the adjoint.  Not available for
        points with a flap."""
        dv = self.dv0 if dv is None else np.asarray(dv, float)
        names = [p.name for p in self.points]
        if points is None:
            pts = list(range(len(self.points)))
        else:
            pts = [names.index(q) if isinstance(q, str) else int(q) for q in np.atleast_1d(points)]
        if method not in ("adjoint", "tangent"):
            raise ValueError("method must be 'adjoint' or 'tangent'")
        n = len(self.param.coordinates(self.shape_dv(dv))[0])
        scal = [f for f in functions if f not in NODAL_FUNCS]
        nodal = [f for f in functions if f in NODAL_FUNCS]
        for f in scal:
            if f not in HYDRO_FUNCS:
                raise ValueError(f"unknown function {f!r}: one of {HYDRO_FUNCS + NODAL_FUNCS}")
        fnames = list(scal) + ([f"{f}:{i}" for f in nodal for i in range(n)] if method == "adjoint" else [])
        req = {ip: set(fnames) for ip in pts}
        for ip in pts:                         # trim angle of the incidence-offset points
            if ip in self.trim_of:
                req.setdefault(self.trim_of[ip], set()).add("alpha")
        x, y, results, vals, grads, defl, _ = self._solve(dv, req, True)
        out = []
        for ip in pts:
            r, p = results[ip], self.points[ip]
            e = dict(name=p.name, mode=p.mode, x=r.x, y=r.y, values={}, dx={}, dy={}, dop={}, ddv={})
            for f in functions:
                if f in NODAL_FUNCS and method == "tangent":
                    continue
                keys = [(ip, f)] if f not in NODAL_FUNCS else [(ip, f"{f}:{i}") for i in range(n)]
                g = [grads[k] for k in keys]
                e["dx"][f] = np.array([gi[0] for gi in g])
                e["dy"][f] = np.array([gi[1] for gi in g])
                e["dop"][f] = np.array([gi[2] for gi in g])
                e["ddv"][f] = np.array([self._chain(grads, defl, k) for k in keys])
                e["values"][f] = np.array([vals[k] for k in keys])
                if f not in NODAL_FUNCS:
                    for key in ("dx", "dy", "dop", "ddv", "values"):
                        e[key][f] = e[key][f][0]
            if nodal and method == "tangent":
                self._tangent_nodal(ip, r, nodal, grads, defl, e)
            out.append(e)
        return out

    def _tangent_nodal(self, ip, r, nodal, grads, defl, e):
        """d(Cp, Cf)/d(design vector) of point ip by tangent passes (one per variable)."""
        from .solver import XFoilTangent
        p = self.points[ip]
        if ip in defl:
            raise NotImplementedError("method='tangent' is not available for points with a flap")
        if getattr(self, "_tangent", None) is None:
            self._tangent = XFoilTangent()
        tg = self._tangent
        op = dict(cl=p.cl) if p.mode == "cl" else dict(alpha=p.alpha)
        rt, _ = solve_point(tg, p, r.x, r.y, r.alpha, self.hydro.tol)
        if rt is None:
            raise EvaluationError(f"tangent library: XFOIL failed at {p.name}")
        Jx, Jy = self.param.jacobian()
        c, dc = self._chord_now
        dtrim = self._chain(grads, defl, (self.trim_of[ip], "alpha")) if ip in self.trim_of else None
        cols = {f: np.zeros((len(r.x), self.ndv)) for f in nodal}
        for k in range(self.n_shape):
            dre = p.re / c * dc[k] if p.speed is not None else 0.0
            dop = dtrim[k] if dtrim is not None else 0.0
            t = tg.tangent(r.x, r.y, dx=Jx[:, k], dy=Jy[:, k], dop=dop, dflow=(dre, 0.0, 0.0),
                           kiter=self.hydro.kiter if "cf" not in nodal else max(self.hydro.kiter, self.hydro.cf_kiter),
                           **op)
            for f in nodal:
                cols[f][:, k] = t[f]
        for f in nodal:
            e["values"][f] = getattr(rt, f).copy()
            e["ddv"][f] = cols[f]
            e["dx"][f] = e["dy"][f] = e["dop"][f] = None

    def write_sensitivities(self, out_dir, dv=None, functions=("cl", "cd", "cp"), points=None, prefix="sens",
                            method="adjoint"):
        """Writes ``sensitivities`` to out_dir, per point: <prefix>_<point>.csv
        (one row per panel node: x, y, Cp, Cf, then dF/dx and dF/dy of the
        scalar functions and the diagonal dCp_i/dx_i, dCp_i/dy_i of the nodal
        ones (adjoint method), and dF_i/d(design variable) of the nodal ones)
        and <prefix>_<point>.npz (all arrays, including the full (n, n)
        matrices of 'cp' / 'cf' and the gradients w.r.t. the design vector);
        and <prefix>_dv.csv: gradients w.r.t. the design variables, one row
        per variable, one column per (function, point) for the scalar functions.
        Returns the list of written paths."""
        import re
        os.makedirs(out_dir, exist_ok=True)
        sens = self.sensitivities(dv, functions, points, method=method)
        dv = self.dv0 if dv is None else np.asarray(dv, float)
        res = self.evaluate(dv, gradients=False)["results"]
        names = [p.name for p in self.points]
        paths, dv_cols = [], {}
        for e in sens:
            ip = names.index(e["name"])
            tag = re.sub(r"[^A-Za-z0-9_.+-]+", "_", e["name"])
            cols = {"x": e["x"], "y": e["y"], "cp": res[ip].cp, "cf": res[ip].cf}
            for f in functions:
                if f in NODAL_FUNCS:
                    if e["dx"][f] is not None:
                        cols[f"d{f}_i/dx_i"] = np.diag(e["dx"][f])
                        cols[f"d{f}_i/dy_i"] = np.diag(e["dy"][f])
                    for k, nm in enumerate(self.names):
                        cols[f"d{f}/d{nm}"] = e["ddv"][f][:, k]
                else:
                    cols[f"d{f}/dx"], cols[f"d{f}/dy"] = e["dx"][f], e["dy"][f]
                    dv_cols[f"{f}@{e['name']}"] = e["ddv"][f]
            path = os.path.join(out_dir, f"{prefix}_{tag}.csv")
            op = "alpha (deg)" if e["mode"] == "alpha" else "CL"
            head = (f"XFOILAD sensitivities, point {e['name']}: derivatives w.r.t. the panel node coordinates "
                    f"at fixed {op}\n" + ",".join(["i"] + list(cols)))
            np.savetxt(path, np.column_stack([np.arange(len(e["x"]))] + list(cols.values())),
                       delimiter=",", header=head, fmt="%.10g")
            paths.append(path)
            arrays = dict(x=e["x"], y=e["y"], cp=res[ip].cp, cf=res[ip].cf, dv=dv)
            for f in functions:
                for key in ("values", "dx", "dy", "dop", "ddv"):
                    if e[key][f] is not None:
                        arrays[f"{key}_{f}"] = e[key][f]
            path = os.path.join(out_dir, f"{prefix}_{tag}.npz")
            np.savez(path, **arrays)
            paths.append(path)
        if dv_cols:
            path = os.path.join(out_dir, f"{prefix}_dv.csv")
            with open(path, "w") as fh:
                fh.write(",".join(["variable", "value"] + list(dv_cols)) + "\n")
                for k, nm in enumerate(self.names):
                    fh.write(",".join([nm, f"{dv[k]:.10g}"] + [f"{c[k]:.10g}" for c in dv_cols.values()]) + "\n")
            paths.append(path)
        return paths

    # --------------------------------------------------------- reporting
    def report(self, dv, title=""):
        ev = self.evaluate(dv, gradients=False)
        d = self.shape_dv(dv)
        flap = self.has_flap or any(c.kind == "cmh" for c in self.constraints)
        cmh = self.hinge_moments(ev) if flap else None
        lines = [title] if title else []
        for p in self.points:
            lines.append(f"{p.name:>14s}: {p.conditions}")
        h = self.hydro
        lines.append(f"{'point':>14s} {'alpha':>7s} {'CL':>7s} {'CD':>9s} {'CDf':>9s} {'CM':>8s} "
                     f"{'L/D':>7s} {'Cpmin':>7s} {'Cfmin':>9s} {'xtr_up':>6s} {'xtr_lo':>6s}"
                     + (f" {'flap':>7s} {'Cmh':>8s}" if flap else ""))
        for ip, p in enumerate(self.points):
            r = ev["results"][ip]
            cfm = float(r.cf[cf_mask(r, h.cfmin_xrange, h.cfmin_side)].min())     # true min, on the cfmin nodes
            lines.append(f"{p.name:>14s} {r.alpha:7.3f} {r.cl:7.4f} {r.cd:9.6f} {r.cdf:9.6f} {r.cm:8.4f} "
                         f"{r.cl / r.cd:7.1f} {r.cpmin:7.3f} {cfm:9.2e} {r.xtr_top:6.3f} {r.xtr_bot:6.3f}"
                         + (f" {ev['flaps'][ip]:7.3f} {cmh[ip]:8.5f}" if flap else ""))
        if flap:
            yh = self.hinge_y(dv)
            lines.append(f"  flap: hinge x/c={self.hinge_x(dv):.4g}{' (design variable)' if self.hinge_var else ''}, "
                         f"y/t={self.flap_hinge[1]:g} "
                         f"(y={yh:.4f}), blend={self.flap_blend:g}")
        props = geomfuncs.section_properties(ev["x"], ev["y"])
        tmax = geomfuncs.max_thickness(self.param, d)[0][0]
        lines.append(f"  geometry: t/c max={tmax:.4f}  area={props['area']:.5f}  Ixx={props['Ixx']:.3e}  "
                     f"Iyy={props['Iyy']:.3e}")
        for c, v in zip(self.constraints, ev["cvals"]):
            lo = "-inf" if c.lower is None else f"{c.lower:g}"
            hi = "+inf" if c.upper is None else f"{c.upper:g}"
            vv = f"{v[0]:.5g}" if len(v) == 1 else f"min={v.min():.4g} max={v.max():.4g}"
            lines.append(f"  constraint {c.name:>22s}: {lo} <= {vv} <= {hi}")
        return "\n".join(lines)

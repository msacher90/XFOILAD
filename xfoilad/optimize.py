"""Gradient-based optimization driver (SciPy SLSQP) for a ``Problem``."""
import json
import time
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize

from .problem import EvaluationError


@dataclass
class OptResult:
    dv0: np.ndarray
    dv: np.ndarray
    success: bool
    message: str
    nit: int
    history: list = field(default_factory=list)
    time: float = 0.0

    def to_json(self, path, problem=None):
        d = dict(dv0=self.dv0.tolist(), dv=self.dv.tolist(), success=bool(self.success),
                 message=self.message, nit=int(self.nit), time=self.time, history=self.history)
        if problem is not None:
            d["names"] = list(getattr(problem, "names", getattr(problem.param, "names", [])))
            d["points"] = [dict(name=p.name, cl=p.cl, alpha=p.alpha, re=p.re, mach=p.mach,
                                ncrit=p.ncrit, tu=p.tu, xtr=list(p.xtr), weight=p.weight,
                                flap=p.flap, flap_var=p.flap_var, alpha_from=p.alpha_from, dalpha=p.dalpha,
                                in_objective=p.in_objective, speed=p.speed)
                           for p in problem.points]
        with open(path, "w") as f:
            json.dump(d, f, indent=1)


def design_scales(problem, dv_scale):
    """Scaling of each design variable (``Problem.dv_scales``: dv_scale for
    the shape, flap_scale for the flap deflections)."""
    if hasattr(problem, "dv_scales"):
        return np.asarray(problem.dv_scales(dv_scale), float)
    return np.full(len(problem.dv0), float(dv_scale))


def point_summary(problem, out):
    """Per-point results stored in the history."""
    flaps = out.get("flaps")
    cmh = problem.hinge_moments(out) if flaps is not None and hasattr(problem, "hinge_moments") else None
    if cmh is None and hasattr(problem, "base") and flaps is not None:     # CVaRProblem
        cmh = problem.base.hinge_moments(dict(out, dv=out["dv"][:problem.base.ndv]))
    pts = {}
    for ip, p in enumerate(problem.points):
        r = out["results"][ip]
        pts[p.name] = dict(cl=r.cl, cd=r.cd, cm=r.cm, alpha=r.alpha, cpmin=r.cpmin)
        if flaps is not None and (p.flap_var or p.flap):
            pts[p.name]["flap"] = float(flaps[ip])
            pts[p.name]["cmh"] = float(cmh[ip])
    return pts


def optimize(problem, maxiter=60, ftol=1e-6, dv_scale=0.01, verbose=True, callback=None, method="slsqp",
             feas_tol=1e-4, **ipopt_kw):
    """Minimizes the problem's objective with SLSQP (default) or IPOPT
    (``method="ipopt"``, see ``xfoilad.ipopt.optimize_ipopt``).

    Design variables are scaled as z = (dv - dv0) / dv_scale (flap
    deflections: z = (delta - delta0) / problem.flap_scale). Hydrodynamic
    functions are evaluated by XFOIL, their gradients by the adjoint.  When
    XFOIL fails to converge at a trial design, a penalized objective is
    returned so that the line search backtracks."""
    if method == "ipopt":
        from .ipopt import optimize_ipopt
        return optimize_ipopt(problem, maxiter=maxiter, tol=ftol, dv_scale=dv_scale, verbose=verbose,
                              callback=callback, **ipopt_kw)
    if method != "slsqp":
        raise ValueError(f"unknown method {method!r}")
    dv0 = np.array(problem.dv0, float)
    scale = design_scales(problem, dv_scale)
    to_dv = lambda z: dv0 + scale * np.asarray(z)
    history = []
    state = dict(last_good=None, last_grad=None, nfev=0, t0=time.time())

    def ev(z, grad=True):
        dv = to_dv(z)
        try:
            out = problem.evaluate(dv, gradients=grad)
            state["last_good"] = out
            if out["cgrads"] is not None:
                state["last_grad"] = out
            return out
        except EvaluationError as e:
            if verbose:
                print(f"    ! {e} -> penalized")
            return None

    # constraint layout: list of (constraint index, side, sign) with g >= 0
    layout = []
    for ic, c in enumerate(problem.constraints):
        if c.lower is not None and c.upper is not None and c.lower == c.upper:
            layout.append((ic, "eq"))
        else:
            if c.lower is not None:
                layout.append((ic, "lo"))
            if c.upper is not None:
                layout.append((ic, "hi"))

    # values without the adjoint (line-search trial points); the gradients of
    # an evaluated design restart from its converged BL states (Problem cache)
    def f(z):
        out = ev(z, grad=False)
        state["nfev"] += 1
        if out is None:
            return 1e3
        return out["f"]

    def df(z):
        out = ev(z, grad=True)
        if out is None or out["grad"] is None:     # last design evaluated with its gradients
            lg = state["last_grad"]
            return np.zeros_like(z) if lg is None or lg["grad"] is None else lg["grad"] * scale
        return out["grad"] * scale

    def make_con(ic, side):
        c = problem.constraints[ic]

        def g(z):
            out = ev(z, grad=False)
            if out is None:
                return -np.ones_like(state["last_good"]["cvals"][ic])
            v = out["cvals"][ic]
            if side == "lo" or side == "eq":
                return (v - c.lower) / c.scale
            return (c.upper - v) / c.scale

        def dg(z):
            out = ev(z, grad=True)
            if out is None:             # last design evaluated with its gradients
                out = state["last_grad"]
            J = out["cgrads"][ic] * scale
            return (J if side in ("lo", "eq") else -J) / c.scale

        return dict(type="eq" if side == "eq" else "ineq", fun=g, jac=dg)

    cons = [make_con(ic, side) for ic, side in layout]
    zlb = (problem.lb - dv0) / scale
    zub = (problem.ub - dv0) / scale

    z0 = np.zeros_like(dv0)
    out0 = problem.evaluate(dv0)
    if verbose:
        print(problem.report(dv0, "=== initial design ==="))
        print(f"\n{'it':>4s} {'objective':>11s} {'max viol':>9s} {'nfev':>5s} {'time':>7s}")

    def max_violation(out):
        mv = 0.0
        for ic, side in layout:
            c = problem.constraints[ic]
            v = out["cvals"][ic]
            if side in ("lo", "eq"):
                mv = max(mv, float(np.max((c.lower - v) / c.scale)))
            if side in ("hi", "eq"):
                mv = max(mv, float(np.max((v - c.upper) / c.scale)))
        return max(mv, 0.0)

    def record(z, it):
        out = ev(z, grad=False)
        if out is None:
            return
        h = dict(it=it, f=float(out["f"]), max_viol=max_violation(out), nfev=state["nfev"],
                 time=time.time() - state["t0"], dv=to_dv(z).tolist(),
                 points=point_summary(problem, out))
        history.append(h)
        if verbose:
            print(f"{it:4d} {h['f']:11.6f} {h['max_viol']:9.2e} {h['nfev']:5d} {h['time']:7.1f}s")
        if callback is not None:
            callback(h)

    record(z0, 0)
    it_counter = [0]

    def cb(z):
        it_counter[0] += 1
        record(z, it_counter[0])

    res = minimize(f, z0, jac=df, method="SLSQP", bounds=list(zip(zlb, zub)), constraints=cons,
                   callback=cb, options=dict(maxiter=maxiter, ftol=ftol, disp=False))
    message = str(res.message)
    if ev(res.x, grad=False) is None and history:  # XFOIL fails at the last step: last converged iterate
        res.x = (np.asarray(history[-1]["dv"]) - dv0) / scale
        message += " (last converged iterate)"
    dv = to_dv(res.x)
    if not history or not np.array_equal(history[-1]["dv"], dv):    # last SLSQP step gets no callback
        record(res.x, it_counter[0] + 1)
    # at the iteration limit SLSQP returns its last iterate, which may be worse
    # or infeasible: return the best feasible iterate instead
    feasible = [h for h in history if h["max_viol"] <= feas_tol]
    if feasible:
        best = min(feasible, key=lambda h: h["f"])
        if history[-1]["max_viol"] > feas_tol or best["f"] < history[-1]["f"]:
            dv = np.asarray(best["dv"])
            message += f" (best feasible iterate, {best['it']})"
    if verbose:
        print(f"\nSLSQP: {res.message} ({res.nit} iterations, {state['nfev']} evaluations, "
              f"{time.time() - state['t0']:.1f} s)")
        print(problem.report(dv, "=== optimized design ==="))
    return OptResult(dv0=dv0, dv=dv, success=bool(res.success), message=message, nit=int(res.nit),
                     history=history, time=time.time() - state["t0"])

"""Interior-point driver (IPOPT, limited-memory quasi-Newton) for a ``Problem``.

Calls the C interface of the official IPOPT binaries through ctypes, so that
no compiled Python extension is needed.  The DLL is looked up in, in order:
the ``dll`` argument, the ``XFOILAD_IPOPT_DLL`` environment variable, the
copy in ``xfoilad/lib/ipopt`` (Windows: DLLs of the official Ipopt 3.14.20 binaries), and
the system search path (``ipopt-3.dll`` / ``libipopt.so`` / ``libipopt.dylib``).

When XFOIL fails to converge at a trial design, the evaluation callback
returns false: IPOPT then cuts the step (no penalized value is needed)."""
import ctypes as C
import os
import sys
import time

import numpy as np

from .optimize import OptResult, design_scales, point_summary
from .problem import EvaluationError

_NUM, _IDX = C.c_double, C.c_int
_P_NUM, _P_IDX = C.POINTER(_NUM), C.POINTER(_IDX)
_EVAL_F = C.CFUNCTYPE(C.c_bool, _IDX, _P_NUM, C.c_bool, _P_NUM, C.c_void_p)
_EVAL_GRAD_F = C.CFUNCTYPE(C.c_bool, _IDX, _P_NUM, C.c_bool, _P_NUM, C.c_void_p)
_EVAL_G = C.CFUNCTYPE(C.c_bool, _IDX, _P_NUM, C.c_bool, _IDX, _P_NUM, C.c_void_p)
_EVAL_JAC_G = C.CFUNCTYPE(C.c_bool, _IDX, _P_NUM, C.c_bool, _IDX, _IDX, _P_IDX, _P_IDX, _P_NUM, C.c_void_p)
_EVAL_H = C.CFUNCTYPE(C.c_bool, _IDX, _P_NUM, C.c_bool, _NUM, _IDX, _P_NUM, C.c_bool, _IDX, _P_IDX, _P_IDX,
                      _P_NUM, C.c_void_p)
_INTERMEDIATE = C.CFUNCTYPE(C.c_bool, _IDX, _IDX, _NUM, _NUM, _NUM, _NUM, _NUM, _NUM, _NUM, _NUM, _IDX,
                            C.c_void_p)

STATUS = {0: "Solve_Succeeded", 1: "Solved_To_Acceptable_Level", 2: "Infeasible_Problem_Detected",
          3: "Search_Direction_Becomes_Too_Small", 4: "Diverging_Iterates", 5: "User_Requested_Stop",
          6: "Feasible_Point_Found", -1: "Maximum_Iterations_Exceeded", -2: "Restoration_Failed",
          -3: "Error_In_Step_Computation", -4: "Maximum_CpuTime_Exceeded", -5: "Maximum_WallTime_Exceeded",
          -10: "Not_Enough_Degrees_Of_Freedom", -11: "Invalid_Problem_Definition", -12: "Invalid_Option",
          -13: "Invalid_Number_Detected", -100: "Unrecoverable_Exception", -101: "NonIpopt_Exception_Thrown",
          -102: "Insufficient_Memory", -199: "Internal_Error"}

_lib = None


def load_library(dll=None):
    global _lib
    if _lib is not None and dll is None:
        return _lib
    path = dll or os.environ.get("XFOILAD_IPOPT_DLL")
    bundled = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib", "ipopt", "ipopt-3.dll")
    if path is None and sys.platform == "win32" and os.path.exists(bundled):
        path = bundled
    if path is None:
        path = {"win32": "ipopt-3.dll", "darwin": "libipopt.dylib"}.get(sys.platform, "libipopt.so")
    elif sys.platform == "win32":                   # dependent DLLs (MUMPS, Intel runtime) sit next to it
        os.add_dll_directory(os.path.dirname(os.path.abspath(path)))
    lib = C.CDLL(path)
    lib.CreateIpoptProblem.restype = C.c_void_p
    lib.CreateIpoptProblem.argtypes = [_IDX, _P_NUM, _P_NUM, _IDX, _P_NUM, _P_NUM, _IDX, _IDX, _IDX,
                                       _EVAL_F, _EVAL_G, _EVAL_GRAD_F, _EVAL_JAC_G, _EVAL_H]
    lib.FreeIpoptProblem.argtypes = [C.c_void_p]
    lib.AddIpoptStrOption.argtypes = [C.c_void_p, C.c_char_p, C.c_char_p]
    lib.AddIpoptNumOption.argtypes = [C.c_void_p, C.c_char_p, _NUM]
    lib.AddIpoptIntOption.argtypes = [C.c_void_p, C.c_char_p, _IDX]
    for f in (lib.AddIpoptStrOption, lib.AddIpoptNumOption, lib.AddIpoptIntOption):
        f.restype = C.c_bool
    lib.SetIntermediateCallback.argtypes = [C.c_void_p, _INTERMEDIATE]
    lib.SetIntermediateCallback.restype = C.c_bool
    lib.GetIpoptCurrentIterate.argtypes = [C.c_void_p, C.c_bool, _IDX, _P_NUM, _P_NUM, _P_NUM, _IDX,
                                           _P_NUM, _P_NUM]
    lib.GetIpoptCurrentIterate.restype = C.c_bool
    lib.IpoptSolve.argtypes = [C.c_void_p, _P_NUM, _P_NUM, _P_NUM, _P_NUM, _P_NUM, _P_NUM, C.c_void_p]
    lib.IpoptSolve.restype = C.c_int
    _lib = lib
    return lib


def _arr(a):
    a = np.ascontiguousarray(a, dtype=np.float64)
    return a, a.ctypes.data_as(_P_NUM)


def optimize_ipopt(problem, maxiter=60, tol=1e-6, dv_scale=0.01, verbose=True, callback=None, dll=None,
                   options=None):
    """Minimizes the problem's objective with IPOPT (L-BFGS Hessian).

    Same scaling as ``optimize`` (z = (dv - dv0) / dv_scale) and same history
    format; two-sided constraints are passed as such.  ``options`` are extra
    IPOPT options, e.g. {"mu_strategy": "monotone"}."""
    lib = load_library(dll)
    dv0 = np.array(problem.dv0, float)
    scale = design_scales(problem, dv_scale)
    n = len(dv0)
    to_dv = lambda z: dv0 + scale * np.asarray(z)
    history = []
    state = dict(nfev=0, nfail=0, t0=time.time())

    def ev(z, grad=True):
        try:
            return problem.evaluate(to_dv(z), gradients=grad)
        except EvaluationError as e:
            state["nfail"] += 1
            if verbose:
                print(f"    ! {e} -> step rejected")
            return None

    # constraint rows: scaled values, l <= g <= u, both sides at once
    gl, gu, rows = [], [], []
    out0 = problem.evaluate(dv0)
    for ic, c in enumerate(problem.constraints):
        for k in range(len(out0["cvals"][ic])):
            rows.append((ic, k))
            gl.append(-2e19 if c.lower is None else c.lower / c.scale)
            gu.append(2e19 if c.upper is None else c.upper / c.scale)
    m = len(rows)

    def gvals(out):
        return np.array([out["cvals"][ic][k] / problem.constraints[ic].scale for ic, k in rows])

    def gjac(out):
        return np.array([out["cgrads"][ic][k] / problem.constraints[ic].scale for ic, k in rows]) * scale

    def zvec(n_, x):
        return np.ctypeslib.as_array(x, shape=(n_,)).copy()

    # values are evaluated without the adjoint (trial points of the line
    # search); the gradients of an evaluated design restart from its
    # converged BL states, kept in the Problem cache (no new XFOIL solve)
    @_EVAL_F
    def eval_f(n_, x, new_x, obj, _):
        state["nfev"] += 1
        out = ev(zvec(n_, x), grad=False)
        if out is None:
            return False
        obj[0] = out["f"]
        return True

    @_EVAL_GRAD_F
    def eval_grad_f(n_, x, new_x, grad, _):
        out = ev(zvec(n_, x))
        if out is None:
            return False
        np.ctypeslib.as_array(grad, shape=(n_,))[:] = out["grad"] * scale
        return True

    @_EVAL_G
    def eval_g(n_, x, new_x, m_, g, _):
        out = ev(zvec(n_, x), grad=False)
        if out is None:
            return False
        np.ctypeslib.as_array(g, shape=(m_,))[:] = gvals(out)
        return True

    @_EVAL_JAC_G
    def eval_jac_g(n_, x, new_x, m_, nele, irow, jcol, values, _):
        if not values:                              # structure: dense, row-major
            ii, jj = np.meshgrid(np.arange(m_), np.arange(n_), indexing="ij")
            np.ctypeslib.as_array(irow, shape=(nele,))[:] = ii.ravel()
            np.ctypeslib.as_array(jcol, shape=(nele,))[:] = jj.ravel()
            return True
        out = ev(zvec(n_, x))
        if out is None:
            return False
        np.ctypeslib.as_array(values, shape=(nele,))[:] = gjac(out).ravel()
        return True

    @_EVAL_H
    def eval_h(*args):                              # never called with the L-BFGS approximation
        return False

    lb, ub = (problem.lb - dv0) / scale, (problem.ub - dv0) / scale
    keep = [_arr(a) for a in (lb, ub, np.array(gl), np.array(gu))]
    handle = lib.CreateIpoptProblem(n, keep[0][1], keep[1][1], m, keep[2][1], keep[3][1], m * n, 0, 0,
                                    eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h)
    if not handle:
        raise RuntimeError("CreateIpoptProblem failed")
    opts = {"hessian_approximation": "limited-memory", "max_iter": maxiter, "tol": tol,
            "acceptable_tol": 10 * tol, "print_level": 5 if verbose else 0, "sb": "yes",
            "nlp_scaling_method": "none", "mu_strategy": "adaptive"}
    opts.update(options or {})
    for k, v in opts.items():
        kb = k.encode()
        ok = (lib.AddIpoptStrOption(handle, kb, v.encode()) if isinstance(v, str) else
              lib.AddIpoptIntOption(handle, kb, v) if isinstance(v, int) else
              lib.AddIpoptNumOption(handle, kb, float(v)))
        if not ok:
            raise ValueError(f"IPOPT option {k}={v!r} rejected")

    def max_violation(out):
        g = gvals(out)
        return float(max(0.0, np.max(np.maximum(np.array(gl) - g, g - np.array(gu)), initial=0.0)))

    def record(z, it):
        out = ev(z, grad=False)
        if out is None:
            return
        h = dict(it=it, f=float(out["f"]), max_viol=max_violation(out), nfev=state["nfev"],
                 time=time.time() - state["t0"], dv=to_dv(z).tolist(),
                 points=point_summary(problem, out))
        history.append(h)
        if callback is not None:
            callback(h)

    zbuf, zbp = _arr(np.zeros(n))

    @_INTERMEDIATE
    def intermediate(alg_mod, it, obj, inf_pr, inf_du, mu, d_norm, reg, a_du, a_pr, ls, _):
        if alg_mod == 0 and lib.GetIpoptCurrentIterate(handle, False, n, zbp, None, None, 0, None, None):
            record(zbuf.copy(), int(it))
        return True

    lib.SetIntermediateCallback(handle, intermediate)
    z, zp = _arr(np.zeros(n))
    obj = _NUM(0.0)
    try:
        status = lib.IpoptSolve(handle, zp, None, C.byref(obj), None, None, None, None)
    finally:
        lib.FreeIpoptProblem(handle)
    msg = STATUS.get(status, str(status))
    if ev(z, grad=False) is None and history:      # XFOIL does not reconverge there: last good iterate
        z = (np.asarray(history[-1]["dv"]) - dv0) / scale
        msg += " (last converged iterate)"
    dv = to_dv(z)
    if not history or not np.array_equal(history[-1]["dv"], dv):
        record(z, history[-1]["it"] + 1 if history else 0)
    nit = history[-1]["it"] if history else 0
    if verbose:
        print(f"\nIPOPT: {msg} ({nit} iterations, {state['nfev']} evaluations, {state['nfail']} rejected, "
              f"{time.time() - state['t0']:.1f} s)")
        print(problem.report(dv, "=== optimized design ==="))
    res = OptResult(dv0=dv0, dv=dv, success=status in (0, 1), message=msg, nit=nit, history=history,
                    time=time.time() - state["t0"])
    res.nfev, res.nfail = state["nfev"], state["nfail"]
    return res

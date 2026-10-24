"""ctypes interface to the adjoint-enabled XFOIL library (xfoilad.dll).

The library holds one global XFOIL state, so a process uses one ``XFoilAD``
instance; the DLL is copied to a temporary file so that several Python
processes (e.g. parallel runs) do not share it.
"""
import ctypes
import os
import shutil
import sys
import tempfile
from ctypes import POINTER, byref, c_double, c_int
from dataclasses import dataclass, field

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if sys.platform.startswith("win"):
    _LIBNAME = "xfoilad.dll"
elif sys.platform == "darwin":
    _LIBNAME = "libxfoilad.dylib"
else:
    _LIBNAME = "libxfoilad.so"

_dp = POINTER(c_double)
_ip = POINTER(c_int)

#: names of the scalar outputs, in the order of the adjoint seed vector
SCALARS = ("cl", "cd", "cdf", "cm", "alpha")


def _arr(a):
    return np.ascontiguousarray(a, dtype=np.float64)


def _ptr(a):
    return a.ctypes.data_as(_dp)


@dataclass
class BLState:
    """Converged boundary-layer state of an analysis (``Result.state``,
    ``XFoilAD.get_state``): theta, delta*, ctau / amplification, Ue and mass
    defect at every BL station, surface velocities, alpha, CL and transition
    intervals.  ``analyze(..., init=state)`` restarts from it, on the same
    section at another operating point (alpha, CL, Re, Ncrit) or on a nearby
    section with the same number of nodes; a warm start converges in a few
    Newton iterations instead of about ten from scratch."""
    r: np.ndarray = field(repr=False)
    i: np.ndarray = field(repr=False)
    #: number of panel nodes of the section it was computed on
    n: int = 0
    alpha: float = float("nan")
    cl: float = float("nan")
    #: flow conditions of that analysis (``XFoilAD.flow``)
    flow: dict = None
    converged: bool = True

    def copy(self):
        return BLState(self.r.copy(), self.i.copy(), self.n, self.alpha, self.cl,
                       dict(self.flow) if self.flow is not None else None, self.converged)


@dataclass
class Result:
    """Converged (or not) XFOIL operating point."""
    cl: float
    cd: float
    cdf: float
    cm: float
    alpha: float
    xtr_top: float
    xtr_bot: float
    cp: np.ndarray
    converged: bool
    x: np.ndarray = field(repr=False, default=None)
    y: np.ndarray = field(repr=False, default=None)
    #: nodal skin friction coefficient tau_w / q_inf, signed along the BL
    #: marching direction from the stagnation point (< 0: separated); 0 if inviscid
    cf: np.ndarray = field(repr=False, default=None)
    #: converged BL state (``BLState``, restart point of ``analyze(init=...)``);
    #: None if the analysis did not converge or is inviscid
    state: BLState = field(repr=False, default=None)

    @property
    def cdp(self):
        return self.cd - self.cdf

    @property
    def cpmin(self):
        return float(self.cp.min())

    @property
    def cfmin(self):
        return float(self.cf.min())


class XFoilAD:
    """Primal and adjoint XFOIL analysis of an airfoil given by its panel nodes.

    Nodes are ordered counterclockwise: upper TE -> LE -> lower TE, and are
    used as XFOIL panel nodes as is (no repaneling inside XFOIL).
    """

    def __init__(self, lib_path=None):
        src = lib_path or os.path.join(_HERE, "lib", _LIBNAME)
        if not os.path.isfile(src):
            raise FileNotFoundError(f"{src} not found: download xfoilad-win64.zip from "
                                    "https://github.com/msacher90/XFOILAD/releases and extract it into xfoilad/lib")
        fd, self._tmp = tempfile.mkstemp(suffix=os.path.splitext(src)[1])
        os.close(fd)
        shutil.copy2(src, self._tmp)
        self._lib = ctypes.CDLL(self._tmp)
        self._lib.xad_init()
        nmax, nblmax = c_int(), c_int()
        self._lib.xad_limits(byref(nmax), byref(nblmax))
        self.nmax, self.nblmax = nmax.value, nblmax.value
        self._has_state = hasattr(self._lib, "xad_state_size")    # absent from tangent libraries before 18/10/2026
        self.flow = None
        self.set_flow()

    def __del__(self):
        lib = getattr(self, "_lib", None)
        if lib is None:
            return
        handle = lib._handle
        del self._lib
        try:
            if sys.platform.startswith("win"):
                ctypes.windll.kernel32.FreeLibrary(ctypes.c_void_p(handle))
            os.remove(self._tmp)
        except OSError:
            pass

    # ------------------------------------------------------------------ setup
    def set_flow(self, re=1e6, mach=0.0, ncrit=9.0, xtr=(1.0, 1.0), verbose=False):
        """Flow conditions (re <= 0: inviscid). ``ncrit`` may be a pair (top, bottom)."""
        nc = (ncrit, ncrit) if np.isscalar(ncrit) else tuple(ncrit)
        self.flow = dict(re=re, mach=mach, ncrit=nc, xtr=tuple(xtr))
        self._lib.xad_set_flow(byref(c_double(re)), byref(c_double(mach)),
                               byref(c_double(nc[0])), byref(c_double(nc[1])),
                               byref(c_double(xtr[0])), byref(c_double(xtr[1])),
                               byref(c_int(int(verbose))))

    def set_tolerance(self, eps):
        """Newton convergence tolerance on the rms BL change (XFOIL: 1e-4)."""
        self._lib.xad_set_tolerance(byref(c_double(eps)))

    def set_stag_frac(self, frac):
        """Minimum distance of the stagnation point from the panel nodes, as a
        fraction of the panel length (default 0.01)."""
        self._lib.xad_set_stag_frac(byref(c_double(frac)))

    def get_state(self, n=0, converged=True):
        """``BLState`` copy of the BL state of the last analysis (n: number of
        nodes of its section, recorded for the check of ``analyze(init=...)``)."""
        nr, ni = c_int(), c_int()
        self._lib.xad_state_size(byref(nr), byref(ni))
        r, i = np.zeros(nr.value), np.zeros(ni.value, dtype=np.int32)
        self._lib.xad_get_state(_ptr(r), i.ctypes.data_as(_ip))
        # the buffer ends with alpha (rad) and CL, see xad_get_state
        return BLState(r, i, int(n), float(np.degrees(r[-2])), float(r[-1]),
                       dict(self.flow) if self.flow is not None else None, bool(converged))

    def set_state(self, state):
        """Restores a BL state (``BLState``, or the pair (r, i) of earlier
        versions); the next ``analyze(..., warm=True)`` restarts from it.
        ``analyze(..., init=state)`` does both."""
        r, i = (state.r, state.i) if isinstance(state, BLState) else state
        r, i = _arr(r), np.ascontiguousarray(i, dtype=np.int32)
        self._lib.xad_set_state(_ptr(r), i.ctypes.data_as(_ip))

    def set_alpha_guess(self, alpha_deg):
        self._lib.xad_set_alpha_guess(byref(c_double(alpha_deg)))

    def _check(self, x, y):
        x, y = _arr(x), _arr(y)
        if x.shape != y.shape or x.ndim != 1:
            raise ValueError("x and y must be 1-D arrays of the same length")
        if len(x) > self.nmax:
            raise ValueError(f"too many nodes ({len(x)} > {self.nmax})")
        if np.any(np.hypot(np.diff(x), np.diff(y)) < 1e-12):
            raise ValueError("duplicated consecutive nodes")
        return x, y

    # ---------------------------------------------------------------- primal
    def analyze(self, x, y, alpha=None, cl=None, niter=100, warm=False, init=None):
        """Solves one operating point at given ``alpha`` (deg) or given ``cl``.

        ``init``: ``BLState`` to restart from (``Result.state`` of an earlier
        analysis, possibly at other flow conditions or on a nearby section
        with the same number of nodes); ``warm=True`` without ``init``
        restarts from the current BL state of the library.  A failed warm
        start is returned as is (``converged=False``): no automatic retry."""
        x, y = self._check(x, y)
        if (alpha is None) == (cl is None):
            raise ValueError("give exactly one of alpha / cl")
        if init is not None:
            if init.n != len(x):
                raise ValueError(f"init state computed on {init.n} nodes, the section has {len(x)}")
            self.set_state(init)
            warm = True
        mode, op = (1, alpha) if alpha is not None else (2, cl)
        n = len(x)
        res = np.zeros(7)
        cp, cf = np.zeros(n), np.zeros(n)
        conv = c_int()
        self._lib.xad_analyze(byref(c_int(n)), _ptr(x), _ptr(y), byref(c_double(op)),
                              byref(c_int(mode)), byref(c_int(niter)), byref(c_int(int(warm))),
                              _ptr(res), _ptr(cp), _ptr(cf), byref(conv))
        ok = (bool(conv.value) and np.all(np.isfinite(res[:5])) and np.all(np.isfinite(cp))
              and np.all(np.isfinite(cf)))
        st = self.get_state(n) if ok and self.flow["re"] > 0 and self._has_state else None
        return Result(*res, cp=cp, converged=ok, x=x, y=y, cf=cf, state=st)

    # --------------------------------------------------------------- adjoint
    def adjoint(self, x, y, alpha=None, cl=None, seed=None, cp_seed=None, cf_seed=None, kiter=3, flow=False):
        """Gradient of  J = sum(seed[k]*out_k) + sum(cp_seed*cp) + sum(cf_seed*cf)
        w.r.t. (x, y, op).

        Must follow a converged ``analyze`` call on the same airfoil and
        operating point: the adjoint replays ``kiter`` Newton iterations from
        that converged state.  ``seed`` is a dict over SCALARS
        (e.g. {"cd": 1.0}).  Returns (dJ/dx, dJ/dy, dJ/dop) where op is alpha
        (deg) in alpha mode and CL in CL mode; with ``flow=True``, also the
        dict {"re": dJ/dRe, "ncrit": (dJ/dNcrit_top, dJ/dNcrit_bottom)} (the
        Reynolds number and Ncrit of ``set_flow``).
        """
        x, y = self._check(x, y)
        mode, op = (1, alpha) if alpha is not None else (2, cl)
        n = len(x)
        s = np.zeros(len(SCALARS))
        for k, v in (seed or {}).items():
            s[SCALARS.index(k)] = v
        cps = np.zeros(n) if cp_seed is None else _arr(cp_seed)
        cfs = np.zeros(n) if cf_seed is None else _arr(cf_seed)
        gx, gy, gop, gflow = np.zeros(n), np.zeros(n), c_double(), np.zeros(3)
        res = np.zeros(len(SCALARS))
        conv = c_int()
        self._lib.xad_adjoint(byref(c_int(n)), _ptr(x), _ptr(y), byref(c_double(op)),
                              byref(c_int(mode)), byref(c_int(kiter)), _ptr(s), _ptr(cps),
                              _ptr(cfs), _ptr(gx), _ptr(gy), byref(gop), _ptr(gflow), _ptr(res), byref(conv))
        if flow:
            return gx, gy, gop.value, dict(re=float(gflow[0]), ncrit=(float(gflow[1]), float(gflow[2])))
        return gx, gy, gop.value


class XFoilTangent(XFoilAD):
    """Primal and tangent (forward) XFOIL, from a separate library
    (xfoilad_tgt.dll in xfoilad/lib).  ``tangent`` gives the
    derivatives of all the outputs, Cp and Cf included, along one direction
    of the inputs: one pass per design variable gives the full Jacobian of
    the Cp distribution, instead of one adjoint pass per node."""

    def __init__(self, lib_path=None):
        name = _LIBNAME.replace("xfoilad", "xfoilad_tgt")
        super().__init__(lib_path or os.path.join(_HERE, "lib", name))

    def tangent(self, x, y, alpha=None, cl=None, dx=None, dy=None, dop=0.0, dflow=(0.0, 0.0, 0.0), kiter=4):
        """Directional derivative of the outputs for the input perturbation
        (dx, dy, dop, dRe, dNcrit_top, dNcrit_bottom), around the converged
        state of the last ``analyze`` on the same airfoil and operating point.
        Returns dict(cl, cd, cdf, cm, alpha, cp, cf) of derivatives."""
        x, y = self._check(x, y)
        mode, op = (1, alpha) if alpha is not None else (2, cl)
        n = len(x)
        xd = np.zeros(n) if dx is None else _arr(dx)
        yd = np.zeros(n) if dy is None else _arr(dy)
        fl = _arr(dflow)
        outd, cpd, cfd, res = np.zeros(len(SCALARS)), np.zeros(n), np.zeros(n), np.zeros(len(SCALARS))
        conv = c_int()
        self._lib.xad_tangent(byref(c_int(n)), _ptr(x), _ptr(y), byref(c_double(op)), byref(c_int(mode)),
                              byref(c_int(kiter)), _ptr(xd), _ptr(yd), byref(c_double(dop)), _ptr(fl),
                              _ptr(outd), _ptr(cpd), _ptr(cfd), _ptr(res), byref(conv))
        out = dict(zip(SCALARS, outd.tolist()))
        out.update(cp=cpd, cf=cfd)
        return out

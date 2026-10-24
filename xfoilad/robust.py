"""Robust design under an uncertain flow parameter: superquantile (CVaR)
objective over a set of weighted samples (operating points).

The CVaR at level alpha of the drag coefficient over the samples k, with
probabilities p_k, is written as in Rockafellar & Uryasev (2000):

    CVaR_a = min_{t, s}  t + sum_k p_k s_k / (1 - a)
             s.t.  s_k >= c_k - t,  s_k >= 0,

where c_k = f_k / f_ref, f_k the objective function of the base problem at
sample k (CD, or chord * CD with objective='drag').  The auxiliary variables
t (the quantile at the optimum) and s_k are appended to the shape variables,
so that the objective and the constraints are smooth; the gradient of each
c_k comes from one adjoint sweep.  With ``groups``, one CVaR per group of
samples (e.g. one per operating condition, the samples being values of the
uncertain incidence) and the objective is their weighted sum.
``CVaRProblem`` exposes the interface used by ``optimize``.

``incidence_samples`` expands operating points into samples of an uncertain
incidence around their trim angle (specified CL) or their incidence.
"""
from dataclasses import replace as _replace
from types import SimpleNamespace

import numpy as np

from .problem import Constraint


def ln_uniform_samples(lo, hi, n):
    """Gauss-Legendre nodes and weights for a log-uniform variable on [lo, hi]."""
    xg, wg = np.polynomial.legendre.leggauss(n)
    u = np.log(lo) + 0.5 * (xg + 1.0) * (np.log(hi) - np.log(lo))
    return np.exp(u), wg / wg.sum()


def uniform_lobatto_samples(lo, hi, n):
    """Gauss-Lobatto nodes and weights for a uniform variable on [lo, hi]: the
    bounds are nodes (e.g. a cavitation constraint at the extreme incidences),
    and the rule is exact for polynomials of degree 2n - 3."""
    if n < 3:
        raise ValueError("Gauss-Lobatto needs n >= 3")
    P = np.polynomial.legendre.Legendre.basis(n - 1)
    xg = np.concatenate([[-1.0], np.sort(P.deriv().roots().real), [1.0]])
    xg = 0.5 * (xg - xg[::-1])                   # exactly symmetric: the middle node is 0 for odd n
    wg = 2.0 / (n * (n - 1) * P(xg) ** 2)
    return lo + 0.5 * (xg + 1.0) * (hi - lo), wg / wg.sum()


def quantile(values, probs, alpha):
    """Quantile (value at risk) of a discrete distribution."""
    i = np.argsort(values)
    v, p = np.asarray(values)[i], np.asarray(probs)[i]
    return float(v[np.searchsorted(np.cumsum(p), alpha - 1e-12)])


def cvar(values, probs, alpha):
    """Superquantile (conditional value at risk) of a discrete distribution."""
    v, p = np.asarray(values, float), np.asarray(probs, float)
    t = quantile(v, p, alpha)
    return float(t + np.sum(p * np.maximum(v - t, 0.0)) / (1.0 - alpha))


class CVaRProblem:
    """CVaR of the objective function f_k of ``base`` (one value per operating
    point in its objective = per sample of the uncertain parameter), with
    probabilities ``probs``.  ``groups``: list of dict(idx=[point indices],
    probs=[...], weight=w) (e.g. from ``incidence_samples``): the objective is
    then sum_g w_g CVaR_a(f_k, k in g).  ``start``: initial shape (default
    base.dv0); f_ref is the probability-weighted f at the start.  The
    geometric and hydrodynamic constraints of ``base`` are kept, also on its
    points with ``in_objective=False``."""

    def __init__(self, base, probs=None, alpha=0.9, t_max=5.0, groups=None, start=None):
        self.base, self.alpha = base, float(alpha)
        obj = base.objective_points
        pos = {ip: k for k, ip in enumerate(obj)}
        if groups is None:
            if probs is None or len(probs) != len(obj):
                raise ValueError("one probability per operating point in the objective")
            groups = [dict(idx=list(obj), probs=list(probs), weight=1.0)]
        self.groups = [g for g in groups if all(i in pos for i in g["idx"])]
        self.rows = [np.array([pos[i] for i in g["idx"]]) for g in self.groups]
        self.pk = [np.asarray(g["probs"], float) / np.sum(g["probs"]) for g in self.groups]
        self.wg = np.array([float(g.get("weight", 1.0)) for g in self.groups])
        self.G = len(self.groups)
        self.N = sum(len(r) for r in self.rows)
        self.p = self.pk[0] if self.G == 1 else None
        self.n = base.ndv
        d0 = base.dv0 if start is None else np.asarray(start, float)
        fk0 = base.evaluate(d0, gradients=False)["fk"]
        self.f_ref = float(sum(w * p @ fk0[r] for w, p, r in zip(self.wg, self.pk, self.rows)))
        self.cd_ref = self.f_ref
        c0 = fk0 / self.f_ref
        t0 = [quantile(c0[r], p, self.alpha) for r, p in zip(self.rows, self.pk)]
        s0 = np.concatenate([np.maximum(c0[r] - t, 0.0) for r, t in zip(self.rows, t0)])
        x0 = np.concatenate([d0, t0, s0])
        self.dv0 = x0
        self.names = (list(base.names) + (["t"] if self.G == 1 else [f"t{g}" for g in range(self.G)])
                      + [f"s{k}" for k in range(self.N)])
        self.param = SimpleNamespace(dv0=x0, ndv=len(x0), base=base.param, names=self.names)
        self.lb = np.concatenate([base.lb, np.zeros(self.G + self.N)])
        self.ub = np.concatenate([base.ub, np.full(self.G + self.N, t_max)])
        self.constraints = list(base.constraints) + [
            Constraint("cvar_tail", lower=0.0, scale=1.0, name="s_k + t - c_k")]
        self.n_base = len(base.constraints)
        self.points = base.points

    def dv_scales(self, dv_scale):
        return np.concatenate([self.base.dv_scales(dv_scale), np.full(self.G + self.N, float(dv_scale))])

    def split(self, x):
        """(design vector of the base problem, t, s); t is a scalar with one group."""
        x = np.asarray(x, float)
        t = x[self.n:self.n + self.G]
        return x[:self.n], (t[0] if self.G == 1 else t), x[self.n + self.G:]

    def evaluate(self, x, gradients=True):
        x = np.asarray(x, float)
        d, t, s = x[:self.n], x[self.n:self.n + self.G], x[self.n + self.G:]
        out = self.base.evaluate(d, gradients=gradients)
        c = out["fk"] / self.f_ref
        k = 1.0 / (1.0 - self.alpha)
        f, tail, grad, Jt, j0 = 0.0, [], np.zeros(len(x)), [], 0
        for g, (r, p) in enumerate(zip(self.rows, self.pk)):
            sg = s[j0:j0 + len(r)]
            f += self.wg[g] * (t[g] + k * p @ sg)
            tail.append(sg + t[g] - c[r])
            grad[self.n + g] = self.wg[g]
            grad[self.n + self.G + j0:self.n + self.G + j0 + len(r)] = self.wg[g] * k * p
            if gradients:
                J = np.zeros((len(r), len(x)))
                J[:, :self.n] = -out["fk_grad"][r] / self.f_ref
                J[:, self.n + g] = 1.0
                J[np.arange(len(r)), self.n + self.G + j0 + np.arange(len(r))] = 1.0
                Jt.append(J)
            j0 += len(r)
        cv = [cvar(c[r], p, self.alpha) for r, p in zip(self.rows, self.pk)]
        mean = float(sum(w * p @ c[r] for w, p, r in zip(self.wg, self.pk, self.rows)))
        res = dict(out, dv=x, f=f, cvals=list(out["cvals"]) + [np.concatenate(tail)], cd=out["fk"], c=c,
                   cvar=float(self.wg @ cv), mean=mean,
                   var=quantile(c[self.rows[0]], self.pk[0], self.alpha) if self.G == 1 else None)
        if gradients:
            pad = lambda J: np.hstack([np.atleast_2d(J), np.zeros((np.atleast_2d(J).shape[0], self.G + self.N))])
            res["grad"] = grad
            res["cgrads"] = [pad(J) for J in out["cgrads"]] + [np.vstack(Jt)]
        else:
            res["grad"] = res["cgrads"] = None
        return res

    def report(self, x, title=""):
        d, t, s = self.split(x)
        ev = self.evaluate(x, gradients=False)
        q = "" if ev["var"] is None else f"  quantile = {ev['var']:.4f}"
        lines = [self.base.report(d, title),
                 f"  CVaR_{self.alpha:g} = {ev['cvar']:.4f}{q}  mean = {ev['mean']:.4f}"
                 f"  (x f_ref = {self.f_ref:.6g})"]
        return "\n".join(lines)


def incidence_samples(points, dalpha=1.5, n=5, fold_zero=False):
    """Expands each operating point into samples of an uncertain incidence,
    uniform on [-dalpha, dalpha] deg, at the nodes of the Gauss-Lobatto rule
    (the bounds are nodes): a specified-CL point gives its trim point (the
    middle node) plus points at alpha_trim + delta (``alpha_from``); an
    alpha-mode point gives points at alpha + delta.  ``fold_zero``: for a
    symmetric section, a point at CL = 0 (or alpha = 0) is computed at
    alpha = 0 and the positive nodes only, CD(-alpha) = CD(alpha) (folded
    probabilities).  The weight of each sample is the weight of its point
    times its probability, so that a Problem(normalize='global') minimizes
    the mean.  Returns (samples, groups) with groups[i] = dict(name, idx,
    dalpha, probs, weight, objective) for point i."""
    da, pr = uniform_lobatto_samples(-dalpha, dalpha, n)
    mid = int(np.argmin(np.abs(da)))
    out, groups = [], []

    def replace(p, **kw):            # copy of a point (tu given: ncrit is recomputed from it)
        return _replace(p, **(dict(ncrit=None) if p.tu is not None else {}), **kw)

    for p in points:
        g = dict(name=p.name, idx=[], dalpha=[], probs=[], weight=p.weight, objective=p.in_objective)

        def add(q, d, pk):
            out.append(q)
            g["idx"].append(len(out) - 1)
            g["dalpha"].append(float(d))
            g["probs"].append(float(pk))

        zero = (p.mode == "cl" and p.cl == 0.0) or (p.mode == "alpha" and p.alpha == 0.0 and p.alpha_from is None)
        if fold_zero and zero:
            for j in range(len(da)):
                if da[j] < -1e-12:
                    continue
                pk = pr[j] if j == mid else 2.0 * pr[j]
                add(replace(p, cl=None, alpha=float(da[j]), weight=p.weight * pk, name=f"{p.name} a{da[j]:+.2f}"),
                    da[j], pk)
        elif p.mode == "cl":
            add(replace(p, weight=p.weight * pr[mid]), 0.0, pr[mid])
            for j in range(len(da)):
                if j != mid:
                    add(replace(p, cl=None, alpha=None, alpha_from=p.name, dalpha=float(da[j]),
                                weight=p.weight * pr[j], name=f"{p.name} da{da[j]:+.2f}"), da[j], pr[j])
        else:
            for j in range(len(da)):
                add(replace(p, alpha=p.alpha + float(da[j]), weight=p.weight * pr[j],
                            name=p.name if j == mid else f"{p.name} da{da[j]:+.2f}"), da[j], pr[j])
        o = np.argsort(g["dalpha"])
        for key in ("idx", "dalpha", "probs"):
            g[key] = [g[key][i] for i in o]
        groups.append(g)
    return out, groups


def ncrit_samples(points, spread=1.0, n=3):
    """Smoothed transition: each point is replaced by n samples of Ncrit,
    uniform on [Ncrit - spread, Ncrit + spread] (Gauss-Legendre nodes), with
    weights = weight of the point x probability.  The free-transition station
    jumps of XFOIL happen at different shapes for different Ncrit, so the
    mean over the samples is a smoother function of the shape (with spread = 1
    and n = 3, the largest second difference of CD along a design variable is
    about 3 times smaller), at n times the cost; the nominal value changes very
    little.  Points given by their Tu are converted to Ncrit first.
    Returns (samples, groups) as ``incidence_samples``."""
    xg, wg = np.polynomial.legendre.leggauss(n)
    wg = wg / wg.sum()
    out, groups = [], []
    for p in points:
        g = dict(name=p.name, idx=[], dncrit=[], probs=[], weight=p.weight, objective=p.in_objective)
        nt, nb = p.ncrit_pair
        for x, w in zip(xg, wg):
            d = float(spread * x)
            q = _replace(p, tu=None, ncrit=(nt + d, nb + d), weight=p.weight * w,
                         name=p.name if abs(d) < 1e-12 else f"{p.name} N{d:+.2f}")
            out.append(q)
            g["idx"].append(len(out) - 1)
            g["dncrit"].append(d)
            g["probs"].append(float(w))
        groups.append(g)
    return out, groups

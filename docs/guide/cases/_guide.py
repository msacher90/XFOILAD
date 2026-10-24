"""Helpers for the guide only: paths, saving results as JSON, summary of an
optimization.  The example scripts import them on their first line; they are
not part of XFOILAD."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(HERE, "..", "out")
FIGS = os.path.join(HERE, "..", "figs")
sys.path.insert(0, ROOT)
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIGS, exist_ok=True)


def _plain(o):
    if hasattr(o, "tolist"):
        return o.tolist()
    if hasattr(o, "__float__"):
        return float(o)
    return str(o)


def save(name, **data):
    """Writes docs/guide/out/<name>.json (read back to write the guide)."""
    with open(os.path.join(OUT, name + ".json"), "w", encoding="utf-8") as f:
        json.dump(data, f, indent=1, default=_plain, ensure_ascii=False)


def fig(name):
    """Path of a figure of the guide."""
    return os.path.join(FIGS, name)


def points_table(pb, dv):
    """Results of each point of a design (no gradient)."""
    ev = pb.evaluate(dv, gradients=False)
    cmh = pb.hinge_moments(ev)
    rows = []
    for ip, p in enumerate(pb.points):
        r = ev["results"][ip]
        rows.append(dict(name=p.name, alpha=r.alpha, cl=r.cl, cd=r.cd, cdf=r.cdf, cm=r.cm, cpmin=r.cpmin,
                         xtr_top=r.xtr_top, xtr_bot=r.xtr_bot, flap=ev["flaps"][ip], cmh=cmh[ip]))
    cons = [dict(name=c.name, lower=c.lower, upper=c.upper, value=v) for c, v in zip(pb.constraints, ev["cvals"])]
    return dict(f=ev["f"], points=rows, constraints=cons)


def summary(pb, res):
    """Initial / optimized design, convergence, time."""
    return dict(initial=points_table(pb, res.dv0), optimized=points_table(pb, res.dv), dv0=res.dv0, dv=res.dv,
                nit=res.nit, time=res.time, message=res.message, ndv=len(res.dv),
                history=[dict(it=h["it"], f=h["f"], max_viol=h["max_viol"]) for h in res.history],
                nfev=res.history[-1]["nfev"] if res.history else 0)

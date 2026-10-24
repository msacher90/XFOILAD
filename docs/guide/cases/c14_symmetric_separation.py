from _guide import save, fig, summary
import numpy as np
import xfoilad as xa
from xfoilad.problem import cf_mask
from xfoilad.plotting import plot_result

# Rudder: symmetric section, minimum CD on course and in a turn; a manoeuvre
# point (CL = 0.9) outside the objective carries the constraint Cf >= 0


def build(cfmin, cfmin_xrange=(0.0, 1.0)):
    param = xa.BSplineAirfoil.from_naca("0012", n_ctrl=10, symmetric=True)      # 8 variables
    points = [xa.OperatingPoint(cl=0.2, re=1e6, tu=0.1, name="course"),
              xa.OperatingPoint(cl=0.4, re=1e6, tu=0.1, name="turn", weight=0.5),
              xa.OperatingPoint(cl=0.9, re=6e5, tu=0.3, name="manoeuvre", in_objective=False)]
    pb = xa.Problem(param, points, objective="cd", cfmin_xrange=cfmin_xrange)
    if cfmin:
        pb.add_constraint("cfmin", point="manoeuvre", lower=0.0)
    pb.add_constraint("alpha", point="manoeuvre", upper=10.0)
    pb.add_constraint("max_thickness", lower=0.12)
    pb.add_constraint("le_radius", lower=0.008)
    return pb


out = {}
for tag, cfmin, xr in (("free", False, (0.0, 1.0)), ("cf_all", True, (0.0, 1.0)), ("cf_x01", True, (0.1, 1.0))):
    print(f"\n===== {tag} =====")
    pb = build(cfmin, xr)
    res = xa.optimize(pb, maxiter=40)
    print(pb.report(res.dv, "=== optimized ==="))
    ev = pb.evaluate(res.dv, gradients=False)
    r = ev["results"][2]
    sep = r.x[r.cf < 0]                       # nodes where the flow is separated (whole section)
    out[tag] = dict(summary(pb, res), cfmin_all=float(r.cf[cf_mask(r)].min()),
                    sep_x=sep.tolist(), x=r.x, cf=r.cf)
    if tag == "cf_x01":
        plot_result(pb, res, fig("c14.png"))
r0 = build(False).evaluate(build(False).dv0, gradients=False)["results"][2]
out["initial"] = dict(cfmin_all=float(r0.cf[cf_mask(r0)].min()), sep_x=r0.x[r0.cf < 0].tolist(), x=r0.x, cf=r0.cf)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig_, ax = plt.subplots(figsize=(8, 3.2))
for tag, lab, st in (("initial", "NACA 0012", ":"), ("free", "unconstrained", "-"),
                     ("cf_all", r"$C_f\geq0$ everywhere", "-"), ("cf_x01", r"$C_f\geq0$ for $x>0.1$", "--")):
    x, cf = np.array(out[tag]["x"]), np.array(out[tag]["cf"])
    ile = int(np.argmin(x))
    ax.plot(x[:ile + 1], 1e3 * cf[:ile + 1], st, lw=1.2, label=lab)    # upper side (suction side)
ax.axhline(0.0, color="k", lw=0.6)
ax.set_xlabel("x/c (upper side)")
ax.set_ylabel(r"$10^3\,C_f$")
ax.set_ylim(-4, 8)
ax.legend(fontsize=8, ncol=2)
ax.set_title("manoeuvre point, CL = 0.9")
fig_.tight_layout()
fig_.savefig(fig("c14_cf.pdf"))
save("c14", **out)

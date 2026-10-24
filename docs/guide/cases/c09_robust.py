from _guide import save, fig
import numpy as np
import xfoilad as xa

# Uncertain Tu, log-uniform on [0.03, 1] %: 5 Gauss-Legendre samples
tus, p = xa.ln_uniform_samples(0.03, 1.0, 5)
param = xa.BSplineAirfoil.from_naca("2412", n_ctrl=10)
pts = [xa.OperatingPoint(cl=0.5, re=1e6, tu=t, name=f"Tu={t:.3g}%") for t in tus]
base = xa.Problem(param, pts, objective="cd")
base.add_constraint("max_thickness", lower=0.12)
base.add_constraint("cpmin", point="all", lower=-1.0)             # at every sample

rob = xa.CVaRProblem(base, p, alpha=0.8)                          # mean of the 20 % worst cases
res = xa.optimize(rob, maxiter=30)
d = res.dv[:base.ndv]                                             # shape variables

cd0 = np.array([r.cd for r in base.evaluate(base.dv0, False)["results"].values()])
cd1 = np.array([r.cd for r in base.evaluate(d, False)["results"].values()])
for t, a, b in zip(tus, cd0, cd1):
    print(f"Tu = {t:.3f} %  CD {a:.5f} -> {b:.5f}")
print("mean", p @ cd0, "->", p @ cd1, "  CVaR_0.8", xa.cvar(cd0, p, 0.8), "->", xa.cvar(cd1, p, 0.8))
xa.write_dat(fig("c09_robust.dat"), *base.coordinates(d), "c09 robust")

save("c09", tus=tus, probs=p, cd0=cd0, cd1=cd1, mean0=p @ cd0, mean1=p @ cd1,
     cvar0=xa.cvar(cd0, p, 0.8), cvar1=xa.cvar(cd1, p, 0.8), nit=res.nit, time=res.time,
     history=[h["f"] for h in res.history])

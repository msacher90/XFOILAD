from _guide import save, fig
import numpy as np
import xfoilad as xa

# Uncertain incidence around the trim (gusts, orbital velocity of the waves):
# alpha = alpha_trim(CL = 0.5) + dalpha, dalpha uniform on [-1.5, 1.5] deg.
# Gauss-Lobatto: the bounds are nodes, the central node is the trim.
das, q = xa.uniform_lobatto_samples(-1.5, 1.5, 5)
param = xa.BSplineAirfoil.from_naca("2412", n_ctrl=10)
pts = [xa.OperatingPoint(cl=0.5, re=1e6, tu=0.07, name="trim")]
pts += [xa.OperatingPoint(alpha_from="trim", dalpha=da, re=1e6, tu=0.07, name=f"{da:+.2f} deg")
        for i, da in enumerate(das) if i != len(das) // 2]
w = np.r_[q[len(das) // 2], np.delete(q, len(das) // 2)]          # probabilities, in the order of the points

base = xa.Problem(param, pts, objective="cd")
cd0 = np.array([r.cd for r in base.evaluate(base.dv0, False)["results"].values()])
for p, wk, c in zip(pts, w, cd0):            # weights p_k CD_k(d0): the objective is the mean of CD
    p.weight = wk * c
pb = xa.Problem(param, pts, objective="cd")
pb.add_constraint("max_thickness", lower=0.12)
for ip, p in enumerate(pts):                 # no cavitation at the extreme incidences
    if p.alpha_from is not None and abs(p.dalpha) == 1.5:
        pb.add_constraint("cpmin", point=ip, lower=-1.0)

res = xa.optimize(pb, maxiter=40)
print(pb.report(res.dv, "optimized"))
ev0, ev1 = pb.evaluate(res.dv0, False), pb.evaluate(res.dv, False)
r0 = [ev0["results"][i] for i in range(len(pts))]
r1 = [ev1["results"][i] for i in range(len(pts))]
xa.write_dat(fig("c12_robust_alpha.dat"), *pb.coordinates(res.dv), "c12 robust in incidence")
save("c12", dalpha=[p.dalpha if p.alpha_from is not None else 0.0 for p in pts], probs=w,
     cd0=[r.cd for r in r0], cd1=[r.cd for r in r1], cp0=[r.cpmin for r in r0], cp1=[r.cpmin for r in r1],
     cl0=[r.cl for r in r0], cl1=[r.cl for r in r1], alpha0=r0[0].alpha, alpha1=r1[0].alpha,
     mean0=float(w @ [r.cd for r in r0]), mean1=float(w @ [r.cd for r in r1]),
     nit=res.nit, time=res.time, history=[h["f"] for h in res.history])

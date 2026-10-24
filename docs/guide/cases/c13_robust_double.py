from _guide import save, fig
import numpy as np
import xfoilad as xa

# Double uncertainty: incidence (uniform on trim +- 1.5 deg, Gauss-Lobatto)
# and turbulence (log-uniform on [0.03, 1] %, Gauss-Legendre), independent.
das, q = xa.uniform_lobatto_samples(-1.5, 1.5, 3)                  # -1.5, 0, +1.5 deg
tus, p = xa.ln_uniform_samples(0.03, 1.0, 3)
param = xa.BSplineAirfoil.from_naca("2412", n_ctrl=10)
pts, w = [], []
for t, pt in zip(tus, p):
    trim = f"Tu={t:.3g}%"                    # one trim per turbulence level: alpha_trim depends on Tu
    pts.append(xa.OperatingPoint(cl=0.5, re=1e6, tu=t, name=trim))
    w.append(pt * q[1])
    for da, qa in ((das[0], q[0]), (das[2], q[2])):
        pts.append(xa.OperatingPoint(alpha_from=trim, dalpha=da, re=1e6, tu=t, name=f"{trim} {da:+g} deg"))
        w.append(pt * qa)
w = np.array(w)

base = xa.Problem(param, pts, objective="cd")
cd0 = np.array([r.cd for r in base.evaluate(base.dv0, False)["results"].values()])
for pk, wk, c in zip(pts, w, cd0):
    pk.weight = wk * c                       # objective: mean of CD over the joint distribution
pb = xa.Problem(param, pts, objective="cd")
pb.add_constraint("max_thickness", lower=0.12)
for ip, pk in enumerate(pts):
    if pk.alpha_from is not None:            # no cavitation at the extreme incidences, at every Tu
        pb.add_constraint("cpmin", point=ip, lower=-1.0)

res = xa.optimize(pb, maxiter=40)
print(pb.report(res.dv, "optimized"))
ev0, ev1 = pb.evaluate(res.dv0, False), pb.evaluate(res.dv, False)
r0 = [ev0["results"][i] for i in range(len(pts))]
r1 = [ev1["results"][i] for i in range(len(pts))]
xa.write_dat(fig("c13_robust_double.dat"), *pb.coordinates(res.dv), "c13 robust in incidence and Tu")
save("c13", names=[pk.name for pk in pts], tu=[pk.tu for pk in pts],
     dalpha=[pk.dalpha if pk.alpha_from is not None else 0.0 for pk in pts], probs=w,
     cd0=[r.cd for r in r0], cd1=[r.cd for r in r1], cp0=[r.cpmin for r in r0], cp1=[r.cpmin for r in r1],
     mean0=float(w @ [r.cd for r in r0]), mean1=float(w @ [r.cd for r in r1]),
     nit=res.nit, time=res.time, npoints=len(pts), history=[h["f"] for h in res.history])

from _guide import save, fig
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import xfoilad as xa

param = xa.BSplineAirfoil.from_naca("2412", n_ctrl=12)
xf = xa.XFoilAD()                      # shared solver (serial computation)
alphas = np.arange(-2.0, 9.0, 1.0)
polars = {}
for delta in (0.0, 4.0, 8.0):
    rows = []
    for a in alphas:
        pb = xa.Problem(param, [xa.OperatingPoint(alpha=a, re=1e6, flap=delta)],
                        flap_hinge=(0.75, 0.5), solver=xf, n_workers=1, verbose=False)
        try:
            ev = pb.evaluate(pb.dv0, gradients=False)
        except xa.EvaluationError:
            continue                   # point not converged: skipped
        r = ev["results"][0]
        rows.append((a, r.cl, r.cd, r.cm, pb.hinge_moments(ev)[0]))
    polars[delta] = np.array(rows)

fig_, ax = plt.subplots(1, 3, figsize=(12, 3.6))
for delta, P in polars.items():
    lab = f"flap {delta:g}°"
    ax[0].plot(P[:, 0], P[:, 1], "o-", ms=3, label=lab)
    ax[1].plot(P[:, 2], P[:, 1], "o-", ms=3, label=lab)
    ax[2].plot(P[:, 0], P[:, 4], "o-", ms=3, label=lab)
ax[0].set(xlabel="α (°)", ylabel="$C_L$")
ax[1].set(xlabel="$C_D$", ylabel="$C_L$")
ax[2].set(xlabel="α (°)", ylabel="$C_{mh}$ (hinge, / $q c^2$)")
for a in ax:
    a.grid(alpha=0.3)
ax[0].legend()
fig_.tight_layout()
fig_.savefig(fig("flap_polars.pdf"))

save("c03", polars={str(k): v for k, v in polars.items()})

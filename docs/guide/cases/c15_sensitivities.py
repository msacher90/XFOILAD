from _guide import save, fig, OUT
import os
import time
import numpy as np
import xfoilad as xa

# Adjoint sensitivities of CL, CD and of the Cp distribution with respect to the
# node coordinates, at imposed alpha, and with respect to the design variables
param = xa.BSplineAirfoil.from_naca("2412", n_ctrl=10)
pb = xa.Problem(param, [xa.OperatingPoint(alpha=4.0, re=1e6, name="a4")], n_workers=1)

t0 = time.perf_counter()
s = pb.sensitivities(functions=("cl", "cd", "cp"))[0]     # one adjoint pass per function and per node
t_sens = time.perf_counter() - t0
x, n = s["x"], len(s["x"])
print(f"{n} nodes: dCL/dy {s['dy']['cl'].shape}, dCp/dy {s['dy']['cp'].shape}, {t_sens:.1f} s")
print(f"dCL/dalpha = {s['dop']['cl']:.4f} /deg   dCD/dalpha = {s['dop']['cd']:.6f} /deg")
print("dCD/d(design variables):", np.round(s["ddv"]["cd"], 4))

# check of one column of dCp/dy with finite differences
xf, j, h = pb.hydro.solver, 40, 1e-6
e = np.zeros(n)
e[j] = h
rp, rm = xf.analyze(x, s["y"] + e, alpha=4.0), xf.analyze(x, s["y"] - e, alpha=4.0)
fd = (rp.cp - rm.cp) / (2 * h)
err = np.abs(s["dy"]["cp"][:, j] - fd).max() / np.abs(fd).max()
print(f"column {j} of dCp/dy: max relative deviation from finite differences {err:.1e}")

paths = pb.write_sensitivities(os.path.join(OUT, "c15_sens"), functions=("cl", "cd", "cp"))
print("\n".join(os.path.relpath(p, OUT) for p in paths))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ile = int(np.argmin(x))
fig_, axs = plt.subplots(1, 3, figsize=(13, 3.6))
for f, ax in zip(("cl", "cd"), axs[:2]):
    g = s["dy"][f]
    ax.plot(x[:ile + 1], g[:ile + 1], label="upper side")
    ax.plot(x[ile:], g[ile:], label="lower side")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("x/c")
    ax.set_title(rf"$\partial C_{{{f[1].upper()}}}/\partial y_j$")
    ax.legend(fontsize=8)
v = np.abs(s["dy"]["cp"]).max() / 20
im = axs[2].imshow(s["dy"]["cp"], cmap="RdBu_r", vmin=-v, vmax=v)
axs[2].set_title(r"$\partial C_{p,i}/\partial y_j$ (node indices)")
axs[2].set_xlabel("j")
axs[2].set_ylabel("i")
fig_.colorbar(im, ax=axs[2], shrink=0.85)
fig_.tight_layout()
fig_.savefig(fig("c15_sens.pdf"))
save("c15", n=n, t_sens=t_sens, dcl_da=s["dop"]["cl"], dcd_da=s["dop"]["cd"], err=err, node=j,
     files=[os.path.relpath(p, OUT) for p in paths], ndv=pb.ndv)

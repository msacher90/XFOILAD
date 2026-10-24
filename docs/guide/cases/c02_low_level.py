from _guide import save
import time
import numpy as np
import xfoilad as xa

x, y = xa.naca4("0012", 80)                    # 159 nodes, counter-clockwise
xf = xa.XFoilAD()                              # one instance per process
xf.set_flow(re=1e6, ncrit=9.0)
xf.set_tolerance(1e-11)

t0 = time.perf_counter()
r = xf.analyze(x, y, cl=0.5)                   # primal, CL mode
t_primal = time.perf_counter() - t0
t0 = time.perf_counter()
gx, gy, gcl = xf.adjoint(x, y, cl=0.5, seed={"cd": 1.0})   # dCD/dx, dCD/dy, dCD/dCL
t_adj = time.perf_counter() - t0

# check on one upper-side node with central finite differences
i, h = 40, 1e-6
e = np.zeros_like(y)
e[i] = h
fd = (xf.analyze(x, y + e, cl=0.5).cd - xf.analyze(x, y - e, cl=0.5).cd) / (2 * h)
print(f"CD = {r.cd:.6f}  alpha = {r.alpha:.3f}  converged = {r.converged}")
print(f"dCD/dy[{i}] adjoint = {gy[i]:.6e}   FD = {fd:.6e}")
print(f"dCD/dCL = {gcl:.5f}   primal {t_primal:.3f} s, adjoint {t_adj:.3f} s")

save("c02", cd=r.cd, alpha=r.alpha, node=i, adj=gy[i], fd=fd, dcd_dcl=gcl,
     t_primal=t_primal, t_adj=t_adj, n=len(x))

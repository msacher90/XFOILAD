from _guide import save
import time
import numpy as np
import xfoilad as xa

param = xa.BSplineAirfoil.from_naca("2412", n_ctrl=12)          # 19 shape variables
pts = [xa.OperatingPoint(cl=0.4, re=1e6, flap=2.0, flap_var=True, name="cruise"),
       xa.OperatingPoint(cl=0.9, re=1e6, name="high CL")]
pb = xa.Problem(param, pts, objective="cd")                       # + 1 flap variable
pb.add_constraint("cpmin", point=1, lower=-2.0)
pb.add_constraint("cmh", point=0, upper=0.01)
pb.hydro.tol = 1e-11                                              # tight Newton for the FD

t0 = time.perf_counter()
ev = pb.evaluate(pb.dv0)                                          # values + adjoint gradients
t_adj = time.perf_counter() - t0
G = np.vstack([ev["grad"]] + [J for J in ev["cgrads"]])           # objective, cpmin, cmh

t0 = time.perf_counter()
FD = np.zeros_like(G)
for k in range(pb.ndv):                                           # central FD, one variable at a time
    h = 1e-6 if k < pb.n_shape else 1e-4                          # flap in degrees
    e = np.zeros(pb.ndv)
    e[k] = h
    fp, fm = pb.evaluate(pb.dv0 + e, False), pb.evaluate(pb.dv0 - e, False)
    FD[:, k] = (np.r_[fp["f"], np.concatenate(fp["cvals"])] - np.r_[fm["f"], np.concatenate(fm["cvals"])]) / (2 * h)
t_fd = time.perf_counter() - t0

err = np.abs(G - FD).max(axis=1) / np.abs(FD).max(axis=1)
for name, e_ in zip(["objective", "cpmin", "cmh"], err):
    print(f"{name:10s} max relative error = {e_:.1e}")
print(f"{pb.ndv} variables: adjoint {t_adj:.2f} s, finite differences {t_fd:.1f} s")

save("c04", err=err, t_adj=t_adj, t_fd=t_fd, ndv=pb.ndv,
     flap_col=dict(adj=G[:, -1], fd=FD[:, -1]))

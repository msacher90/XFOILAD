from _guide import save, fig, summary
import xfoilad as xa
from xfoilad.plotting import plot_result

# 95 shape variables: same gradient cost as with 15
param = xa.BSplineAirfoil.from_naca("2412", n_ctrl=50)
pb = xa.Problem(param, [xa.OperatingPoint(cl=0.5, re=1e6, tu=0.1)], objective="cd", dv_bounds=0.02)
pb.add_constraint("max_thickness", lower=0.12)
pb.add_constraint("cpmin", point=0, lower=-1.0)
pb.add_constraint("smoothness", lower=-1e-3, upper=1e-3)          # no waviness

res = xa.optimize(pb, maxiter=60, dv_scale=0.002)
plot_result(pb, res, fig("c08.png"))
save("c08", nodes=param.n_nodes, **summary(pb, res))

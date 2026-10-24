from _guide import save, fig, summary
import xfoilad as xa
from xfoilad.plotting import plot_result

# One point: minimum CD at CL = 0.5 with a cavitation constraint (Cpmin >= -0.7)
param = xa.BSplineAirfoil.from_naca("2412", n_ctrl=16)            # 27 variables
pb = xa.Problem(param, [xa.OperatingPoint(cl=0.5, re=1e6)], objective="cd")
pb.add_constraint("max_thickness", lower=0.10)
pb.add_constraint("cpmin", point=0, lower=-0.7)

res = xa.optimize(pb, maxiter=45)
xa.write_dat(fig("c05_optimized.dat"), *pb.coordinates(res.dv), "c05 optimized")
plot_result(pb, res, fig("c05.png"))
save("c05", **summary(pb, res))

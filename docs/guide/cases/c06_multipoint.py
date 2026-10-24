from _guide import save, fig, summary
import xfoilad as xa
from xfoilad.plotting import plot_result

# Three weighted points, hydrodynamic and structural constraints
param = xa.BSplineAirfoil.from_naca("2412", n_ctrl=10)            # 15 variables
points = [xa.OperatingPoint(cl=0.3, re=8e5, weight=1.0, name="cruise"),
          xa.OperatingPoint(cl=0.6, re=8e5, weight=1.0, name="climb"),
          xa.OperatingPoint(cl=0.9, re=6e5, weight=0.5, name="high CL")]
pb = xa.Problem(param, points, objective="cd", dv_bounds=0.04)
I0 = xa.geomfuncs.inertia(param, param.dv0)[0][0]                  # initial inertia
pb.add_constraint("cpmin", point=2, lower=-1.2)                    # cavitation margin
pb.add_constraint("alpha", point=2, upper=7.0)
pb.add_constraint("box", x_range=(0.15, 0.55), y_range=(-0.035, 0.045))   # spar box
pb.add_constraint("inertia", axis="xx", lower=I0)
pb.add_constraint("area", lower=0.075)
pb.add_constraint("le_radius", lower=0.006)

res = xa.optimize(pb, maxiter=50)
plot_result(pb, res, fig("c06.png"))
save("c06", I0=I0, **summary(pb, res))

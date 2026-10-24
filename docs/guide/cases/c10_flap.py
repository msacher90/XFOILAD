from _guide import save, fig, summary
import xfoilad as xa
from xfoilad.plotting import plot_result

# Shape and flap deflections optimized together; at imposed CL, XFOIL solves for alpha
param = xa.BSplineAirfoil.from_naca("2412", n_ctrl=12)
points = [xa.OperatingPoint(cl=0.3, re=1e6, flap=0.0, flap_var=True, name="cruise"),
          xa.OperatingPoint(cl=1.0, re=1e6, flap=0.0, flap_var=True, name="climb")]
pb = xa.Problem(param, points, objective="cd",
                flap_hinge=(0.75, 0.5), flap_bounds=(-5.0, 15.0))
pb.add_constraint("max_thickness", lower=0.12)
pb.add_constraint("cm", point=0, lower=-0.08)
pb.add_constraint("cpmin", point=1, lower=-1.6)

res = xa.optimize(pb, maxiter=40)
print("deflections:", pb.flap_angles(res.dv))
plot_result(pb, res, fig("c10.png"))
save("c10", **summary(pb, res))

from _guide import save
import xfoilad as xa

# B-spline section fitted to a NACA 2412 (12 control points per side)
param = xa.BSplineAirfoil.from_naca("2412", n_ctrl=12)

# Three operating points: imposed CL, imposed alpha, deflected flap
points = [xa.OperatingPoint(cl=0.5, re=1e6, name="CL 0.5"),
          xa.OperatingPoint(alpha=4.0, re=1e6, tu=0.5, name="alpha 4"),
          xa.OperatingPoint(cl=0.5, re=1e6, flap=5.0, name="flap 5 deg")]
pb = xa.Problem(param, points)

print(pb.report(pb.dv0))                       # table of results
ev = pb.evaluate(pb.dv0, gradients=False)      # direct access to the values
r = ev["results"][0]
print(f"alpha = {r.alpha:.3f} deg, CD = {r.cd:.6f}, upper-side transition x/c = {r.xtr_top:.3f}")

save("c01", report=pb.report(pb.dv0), alpha=r.alpha, cd=r.cd, xtr=r.xtr_top)

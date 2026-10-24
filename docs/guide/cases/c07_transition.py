from _guide import save
import xfoilad as xa

param = xa.BSplineAirfoil.from_naca("0012", n_ctrl=12)
xf = xa.XFoilAD()
rows = []
for tu in (0.03, 0.07, 0.2, 0.5, 1.0):
    p = xa.OperatingPoint(cl=0.4, re=1e6, tu=tu)
    pb = xa.Problem(param, [p], solver=xf, n_workers=1, verbose=False)
    r = pb.evaluate(pb.dv0, gradients=False)["results"][0]
    rows.append(dict(tu=tu, ncrit=p.ncrit, xtr_top=r.xtr_top,
                     xtr_bot=r.xtr_bot, cd=r.cd, alpha=r.alpha))
    print(f"Tu={tu:4.2f}%  Ncrit={p.ncrit:5.2f}  "
          f"xtr={r.xtr_top:.3f}/{r.xtr_bot:.3f}  CD={r.cd:.5f}")
save("c07", rows=rows)

"""Result figures: shapes, Cp distributions, convergence history."""
import numpy as np


def plot_result(problem, res, path=None, show=False):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    param = problem.param
    ev0 = problem.evaluate(res.dv0, gradients=False)
    r0 = ev0["results"]
    x0, y0 = ev0["x"], ev0["y"]
    ev1 = problem.evaluate(res.dv, gradients=False)
    r1 = ev1["results"]
    x1, y1 = ev1["x"], ev1["y"]

    npts = len(problem.points)
    fig = plt.figure(figsize=(13, 4.2 + 3.2 * ((npts + 1) // 2)))
    gs = fig.add_gridspec(1 + (npts + 1) // 2, 2)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(x0, y0, "--", color="0.5", lw=1, label="initial")
    ax.plot(x1, y1, "-", color="C0", lw=1.6, label="optimized")
    flaps = ev1.get("flaps")
    for ip, pt in enumerate(problem.points):          # sections seen by the flapped points
        if flaps is not None and (pt.flap_var or pt.flap):
            ax.plot(r1[ip].x, r1[ip].y, "-", lw=0.7, alpha=0.8, label=f"{pt.name}: flap {flaps[ip]:.2f} deg")
    if hasattr(param, "control_polygon"):
        (xc, cu), (_, cl) = param.control_polygon(problem.shape_dv(res.dv) if hasattr(problem, "shape_dv")
                                                  else res.dv)
        ax.plot(xc, cu, "o:", color="C1", ms=3, lw=0.8, label="control points")
        ax.plot(xc, cl, "o:", color="C1", ms=3, lw=0.8)
    for c in problem.constraints:
        if c.kind == "box":
            (xa, xb), (ya, yb) = c.options["x_range"], c.options["y_range"]
            ax.add_patch(plt.Rectangle((xa, ya), xb - xa, yb - ya, fill=False, ec="C3", lw=1.2, label="box"))
    ax.set_aspect("equal")
    ax.set_xlabel("x/c")
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.5), ncol=3, frameon=False)
    ax.set_title("Section")

    ax = fig.add_subplot(gs[0, 1])
    it = [h["it"] for h in res.history]
    ax.plot(it, [h["f"] for h in res.history], "o-", ms=3, color="C0", label="objective (normalized)")
    ax.set_xlabel("iteration")
    ax.set_ylabel("objective")
    ax2 = ax.twinx()
    ax2.semilogy(it, [max(h["max_viol"], 1e-12) for h in res.history], "s--", ms=2, color="C3", lw=0.8,
                 label="max violation")
    ax2.set_ylabel("max constraint violation")
    ax.set_title("Convergence")
    ax.grid(alpha=0.3)

    for ip, p in enumerate(problem.points):
        ax = fig.add_subplot(gs[1 + ip // 2, ip % 2])
        ax.plot(r0[ip].x, r0[ip].cp, "--", color="0.5", lw=1,
                label=f"initial  CD={r0[ip].cd:.5f}  a={r0[ip].alpha:.2f}")
        ax.plot(r1[ip].x, r1[ip].cp, "-", color="C0", lw=1.4,
                label=f"optimized CD={r1[ip].cd:.5f}  a={r1[ip].alpha:.2f}")
        for c in problem.constraints:
            if c.kind == "cpmin" and c.point == ip and c.lower is not None:
                ax.axhline(c.lower, color="C3", lw=0.8, ls=":", label=f"Cpmin limit {c.lower:g}")
        ax.invert_yaxis()
        ax.set_xlabel("x/c")
        ax.set_ylabel("Cp")
        ax.set_title(f"{p.name}   {p.conditions}", fontsize=9)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=130)
    if show:
        plt.show()
    plt.close(fig)

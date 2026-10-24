"""Command line: python -m xfoilad config.json

Runs an optimization described in a JSON file (see examples/*.json) and
writes report.txt, history.json, initial.dat, optimized.dat, result.png in
the output directory.

Optional entries:
  "chord": 0.3 | {"inertia": I, "axis": "xx"} | {"modulus": W}   dimensional chord (m, sizing)
  "nu": 1.19e-6                     kinematic viscosity (m^2/s) for the points given by
                                    "speed" (m/s) or "speed_kn" (knots): Re = V c / nu
  "objective": "drag"               chord * CD;  "normalize": "global" | "point"
  "robust_incidence": {"dalpha": 1.5, "n": 5, "fold_zero": true, "measure": "mean" | "cvar",
                       "level": 0.9}
                                    each point -> Gauss-Lobatto samples of the incidence
                                    (mean or CVaR per point); a constraint on a point
                                    applies to all its samples
  "sensitivities": {"functions": ["cl", "cd", "cp"], "designs": ["initial", "optimized"],
                    "method": "adjoint" | "tangent"}
                                    sensitivities w.r.t. the panel nodes and the design
                                    variables (sens_*.csv / .npz, Problem.write_sensitivities)
                                    in <output>/sensitivities_<design>."""
import json
import os
import sys

import numpy as np

from . import (BSplineAirfoil, CSTAirfoil, ChordFromInertia, ChordFromModulus, CVaRProblem, KNOT, OperatingPoint,
               Problem, incidence_samples, optimize, write_dat)
from .plotting import plot_result


def _chord(spec):
    if spec is None or np.isscalar(spec):
        return spec
    if "inertia" in spec:
        return ChordFromInertia(spec["inertia"], spec.get("axis", "xx"))
    if "modulus" in spec:
        return ChordFromModulus(spec["modulus"])
    return float(spec["fixed"])


def build_problem(cfg, base_dir="."):
    """(problem to optimize, base Problem): they differ with a CVaR objective."""
    a = cfg.get("airfoil", {"naca": "0012"})
    pz = dict(cfg.get("parametrization", {}))
    kind = pz.pop("type", "bspline")
    cls = BSplineAirfoil if kind == "bspline" else CSTAirfoil
    if "dat" in a:
        param = cls.from_dat(os.path.join(base_dir, a["dat"]), **pz)
    else:
        param = cls.from_naca(str(a.get("naca", "0012")), **pz)
    points = []
    for p in cfg["points"]:
        p = dict(p)
        if "speed_kn" in p:
            p["speed"] = p.pop("speed_kn") * KNOT
        points.append(OperatingPoint(**p))
    rob = cfg.get("robust_incidence")
    groups = None
    if rob:
        points, groups = incidence_samples(points, dalpha=rob.get("dalpha", 1.5), n=rob.get("n", 5),
                                           fold_zero=rob.get("fold_zero", False))
    o = cfg.get("optimizer", {})
    fl = cfg.get("flap", {})          # {"hinge": [x/c, y/t], "bounds": [lo, hi], "blend": .., "scale": ..,
    #                                    "hinge_var": true, "hinge_bounds": [lo, hi], "hinge_scale": ..}
    pb = Problem(param, points, objective=cfg.get("objective", "cd"), dv_bounds=o.get("dv_bounds", 0.04),
                 kiter=o.get("kiter", 4), tol=o.get("tol", 1e-5), n_workers=o.get("n_workers"),
                 cpmin_rho=o.get("cpmin_rho", 80.0), cpmin_xrange=tuple(o.get("cpmin_xrange", (0.0, 0.98))),
                 flap_hinge=tuple(fl.get("hinge", (0.75, 0.5))), flap_bounds=tuple(fl.get("bounds", (-10.0, 15.0))),
                 flap_blend=fl.get("blend", 0.03), flap_scale=fl.get("scale", 1.0),
                 hinge_var=bool(fl.get("hinge_var", False)), hinge_bounds=tuple(fl.get("hinge_bounds", (0.55, 0.85))),
                 hinge_scale=fl.get("hinge_scale", 0.01),
                 cfmin_rho=o.get("cfmin_rho", 5e4), cfmin_xrange=tuple(o.get("cfmin_xrange", (0.0, 1.0))),
                 cfmin_side=o.get("cfmin_side", "both"), cf_kiter=o.get("cf_kiter", 6),
                 chord=_chord(cfg.get("chord")), nu=cfg.get("nu"),
                 normalize=cfg.get("normalize", "global" if rob else "point"),
                 warm_start=bool(o.get("warm_start", False)))
    for c in cfg.get("constraints", []):
        c = dict(c)
        kind_c = c.pop("kind")
        pt = c.get("point")
        if groups is not None and pt not in (None, "all"):          # every sample of that point
            g = groups[[gr["name"] for gr in groups].index(pt)] if isinstance(pt, str) else groups[int(pt)]
            for ip in g["idx"]:
                pb.add_constraint(kind_c, **dict(c, point=ip))
        else:
            pb.add_constraint(kind_c, **c)
    if rob and rob.get("measure", "mean") == "cvar":
        return CVaRProblem(pb, alpha=rob.get("level", 0.9), groups=[g for g in groups if g["objective"]]), pb
    return pb, pb


def main(argv=None):
    argv = argv or sys.argv[1:]
    if not argv:
        print(__doc__)
        return 1
    path = argv[0]
    with open(path) as f:
        cfg = json.load(f)
    base = os.path.dirname(os.path.abspath(path))
    out = os.path.join(base, cfg.get("output", "results"))
    os.makedirs(out, exist_ok=True)
    prob, pb = build_problem(cfg, base)
    o = cfg.get("optimizer", {})
    res = optimize(prob, maxiter=o.get("maxiter", 60), ftol=o.get("ftol", 1e-6), dv_scale=o.get("dv_scale", 0.01))
    res.to_json(os.path.join(out, "history.json"), prob)
    d0, d1 = res.dv0[:pb.ndv], res.dv[:pb.ndv]
    with open(os.path.join(out, "report.txt"), "w") as f:
        f.write(prob.report(res.dv0, "=== initial design ===") + "\n\n")
        f.write(prob.report(res.dv, "=== optimized design ===") + "\n")
        if pb.chord_model is not None:
            f.write(f"\nchord: {pb.chord_model} -> {pb.chord(d0)[0]:.5f} m -> {pb.chord(d1)[0]:.5f} m\n")
        f.write(f"\nSLSQP: {res.message}, {res.nit} iterations, {res.time:.1f} s\n")
    write_dat(os.path.join(out, "initial.dat"), *pb.coordinates(d0), "XFOILAD initial")
    write_dat(os.path.join(out, "optimized.dat"), *pb.coordinates(d1), "XFOILAD optimized")
    for ip, p in enumerate(pb.points):                # deflected sections, one per flapped point
        if p.flap_var or p.flap:
            write_dat(os.path.join(out, f"optimized_point{ip}.dat"), *pb.coordinates(d1, point=ip),
                      f"XFOILAD optimized, {p.name}, flap {pb.flap_angles(d1)[ip]:.3f} deg")
    plot_result(prob, res, os.path.join(out, "result.png"))
    sens = cfg.get("sensitivities")
    if sens:
        if isinstance(sens, list):
            sens = {"functions": sens}
        for design in sens.get("designs", ["optimized"]):
            dv = {"initial": d0, "optimized": d1}[design]
            pb.write_sensitivities(os.path.join(out, f"sensitivities_{design}"), dv,
                                   functions=tuple(sens.get("functions", ("cl", "cd", "cp"))),
                                   points=sens.get("points"), method=sens.get("method", "adjoint"))
    print(f"results written in {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

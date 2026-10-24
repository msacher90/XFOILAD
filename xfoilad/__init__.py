"""XFOILAD: adjoint-based airfoil shape optimization with an XFOIL core
differentiated by Tapenade (reverse mode)."""
from .solver import XFoilAD, XFoilTangent, Result, BLState
from .geometry import BSplineAirfoil, CSTAirfoil, FlapDeflection, naca4, read_dat, write_dat, repanel
from .problem import (OperatingPoint, Problem, HydroEvaluator, EvaluationError,
                      ncrit_from_tu, tu_from_ncrit, cavitation_number)
from .optimize import optimize, OptResult
from .robust import (CVaRProblem, ln_uniform_samples, uniform_lobatto_samples, quantile, cvar, incidence_samples,
                     ncrit_samples)
from .sizing import FixedChord, ChordFromInertia, ChordFromModulus, KNOT
from . import geomfuncs

__version__ = "0.1.0"

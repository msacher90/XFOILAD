"""Dimensional chord of the section, as a function of the design variables.

A chord model maps the shape variables ``d`` of a parametrization (unit
chord) to the dimensional chord c (m) and its gradient dc/dd.  With a
``Problem(chord=...)``, the Reynolds number of the operating points given by
their speed is Re = V c / nu, and the objective ``"drag"`` is c * CD (drag
per unit span and dynamic pressure, m).

* ``FixedChord(c)``: constant chord.
* ``ChordFromInertia(I, axis="xx")``: the chord that gives the section the
  second moment of area I (m^4): c = (I / I_unit(d))^(1/4), with I_unit the
  inertia of the unit-chord section.  A thinner section becomes longer.
* ``ChordFromModulus(W)``: same with the section modulus (bending stress):
  c = (W / W_unit(d))^(1/3).
"""
import numpy as np

from . import geomfuncs

#: one knot in m/s
KNOT = 1852.0 / 3600.0


class FixedChord:
    def __init__(self, chord):
        self.chord = float(chord)

    def __call__(self, param, d):
        return self.chord, np.zeros(len(d))

    def __repr__(self):
        return f"FixedChord({self.chord:g} m)"


class ChordFromInertia:
    def __init__(self, inertia, axis="xx"):
        self.inertia, self.axis = float(inertia), axis

    def __call__(self, param, d):
        I, dI = geomfuncs.inertia(param, d, axis=self.axis)
        c = (self.inertia / I[0]) ** 0.25
        return c, -0.25 * c / I[0] * dI[0]

    def __repr__(self):
        return f"ChordFromInertia(I{self.axis} = {self.inertia:g} m^4)"


class ChordFromModulus:
    def __init__(self, modulus):
        self.modulus = float(modulus)

    def __call__(self, param, d):
        W, dW = geomfuncs.section_modulus(param, d)
        c = (self.modulus / W[0]) ** (1.0 / 3.0)
        return c, -c / (3.0 * W[0]) * dW[0]

    def __repr__(self):
        return f"ChordFromModulus(W = {self.modulus:g} m^3)"

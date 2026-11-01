"""Shared helpers for the binding tests."""
import numpy as np

from pyrol import Objective
from pyrol.pyrol.Teuchos import ParameterList


def parameters(step="Trust Region"):
    """A quiet parameter list for a solver that uses the given step type."""
    p = ParameterList()
    p["General"] = ParameterList()
    p["General"]["Output Level"] = 0
    p["Step"] = ParameterList()
    p["Step"]["Type"] = step
    return p


class Quadratic(Objective):
    """f(x) = ||x - center||^2 / 2."""

    def __init__(self, center):
        self.center = np.asarray(center, dtype=float)
        super().__init__()

    def value(self, x, tol):
        return 0.5 * float(np.sum((x.array - self.center) ** 2))

    def gradient(self, g, x, tol):
        g.array[:] = x.array - self.center

    def hessVec(self, hv, v, x, tol):
        hv.array[:] = v.array

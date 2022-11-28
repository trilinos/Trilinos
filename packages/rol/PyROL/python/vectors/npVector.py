import numpy as np
from numpy import linalg as LA
from PyROL.PyROL import ROL


class npVector(ROL.Vector_double_t):
    def __init__(self, dimension=1, default_value=0., values=None):
        if values is None:
            self.values = default_value*np.ones((dimension,))
        else:
            self.values = values
        super().__init__()

    def plus(self, b):
        self.values += b.values

    def scale(self, scale_factor):
        self.values *= scale_factor

    def dot(self, b):
        return np.dot(self.values, b.values)

    def norm(self):
        return LA.norm(self.values)

    def clone(self):
        return npVector(dimension=len(self.values))

    def axpy(self, scale_factor, x):
        ax = x.clone()
        ax.plus(x)
        ax.scale(scale_factor)
        self.plus(ax)

    def dimension(self):
        return len(self.values)

    def setScalar(self, new_value):
        self.values[:] = new_value

    def __getitem__(self, index):
        return self.values[index]

    def __setitem__(self, index, val):
        self.values[index] = val


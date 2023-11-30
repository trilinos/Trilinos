import numpy as np
from numpy import linalg as LA
from pyrol.pyrol import ROL
from pyrol.getTypeName import *


class NumPyVector(getTypeName('Vector')):
    def __init__(self, array=None):
        assert isinstance(array, np.ndarray)
        assert array.ndim == 1
        self.array = array
        super().__init__()

    @staticmethod
    def full(dimension=1, default_value=0.):
        array = np.full((dimension,), fill_value=default_value)
        return NumPyVector(array)

    def plus(self, b):
        self.array += b.array

    def scale(self, scale_factor):
        self.array *= scale_factor

    def dot(self, b):
        return np.vdot(self.array, b.array)

    def norm(self):
        return np.sqrt(self.dot(self))

    def zero(self):
        self.setScalar(0.)

    def clone(self):
        tmp = type(self)(np.full(self.array.shape, fill_value=0.))
        return tmp

    def axpy(self, scale_factor, x):
        ax = x.clone()
        ax.zero()
        ax.plus(x)
        ax.scale(scale_factor)
        self.plus(ax)

    def dimension(self):
        return len(self.array)

    def setScalar(self, new_value):
        self.array[:] = new_value

    def reduce(self, op):
        reductionType = op.reductionType()
        if reductionType == ROL.Elementwise.REDUCE_MIN:
            return self.array.min()
        elif reductionType == ROL.Elementwise.REDUCE_MAX:
            return self.array.max()
        elif reductionType == ROL.Elementwise.REDUCE_SUM:
            return self.array.sum()
        elif reductionType == ROL.Elementwise.REDUCE_AND:
            return np.logical_and.reduce(self.array)
        elif reductionType == ROL.Elementwise.REDUCE_BOR:
            return np.bitwise_or.reduce(self.array)
        else:
            raise NotImplementedError(reductionType)

    def applyUnary(self, op):
        for i in range(self.dimension()):
            self.array[i] = op.apply(self.array[i])

    def applyBinary(self, op, other):
        assert self.dimension() == other.dimension()
        for i in range(self.dimension()):
            self.array[i] = op.apply(self.array[i], other.array[i])

    def __getitem__(self, index):
        return self.array[index]

    def __setitem__(self, index, val):
        self.array[index] = val

    def basis(self, i):
        b = self.clone()
        b.array[:] = 0.
        b.array[i] = 1.
        self._basis = b
        return b

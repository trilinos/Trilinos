import numpy as np
from numpy import linalg as LA
from PyROL.PyROL import ROL
from PyROL.getTypeName import *


class npVector(getTypeName('Vector')):
    def __init__(self, values=None):
        assert isinstance(values, np.ndarray)
        assert values.ndim == 1
        self.values = values
        self.copies = []
        super().__init__()

    def __del__(self):
        for copy in self.copies:
            del copy

    @staticmethod
    def full(dimension=1, default_value=0.):
        values = np.full((dimension,), fill_value=default_value)
        return npVector(values)

    def plus(self, b):
        self.values += b.values

    def scale(self, scale_factor):
        self.values *= scale_factor

    def dot(self, b):
        return np.vdot(self.values, b.values)

    def norm(self):
        return LA.norm(self.values)

    def clone(self):
        tmp = npVector(np.empty_like(self.values))
        self.copies.append(tmp)
        return tmp

    def axpy(self, scale_factor, x):
        ax = x.clone()
        ax.zero()
        ax.plus(x)
        ax.scale(scale_factor)
        self.plus(ax)

    def dimension(self):
        return len(self.values)

    def setScalar(self, new_value):
        self.values[:] = new_value

    def reduce(self, op):
        reductionType = op.reductionType()
        if reductionType == ROL.Elementwise.REDUCE_MIN:
            return self.values.min()
        elif reductionType == ROL.Elementwise.REDUCE_MAX:
            return self.values.max()
        elif reductionType == ROL.Elementwise.REDUCE_SUM:
            return self.values.sum()
        elif reductionType == ROL.Elementwise.REDUCE_AND:
            return np.logical_and.reduce(self.values)
        elif reductionType == ROL.Elementwise.REDUCE_BOR:
            return np.bitwise_or.reduce(self.values)
        else:
            raise NotImplementedError(reductionType)

    def applyUnary(self, op):
        for i in range(self.dimension()):
            self.values[i] = op.apply(self.values[i])

    def applyBinary(self, op, other):
        assert self.dimension() == other.dimension()
        for i in range(self.dimension()):
            self.values[i] = op.apply(self.values[i], other.values[i])

    def __getitem__(self, index):
        return self.values[index]

    def __setitem__(self, index, val):
        self.values[index] = val

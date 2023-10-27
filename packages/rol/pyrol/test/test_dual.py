import numpy as np
from numpy import linalg as LA
from PyROL.getTypeName import getTypeName
from PyROL import Vector, Reduced_Objective_SimOpt, Objective_SimOpt, Constraint_SimOpt
from PyROL.vectors import npVector

# Testing 5 implementations of ROL::Vector:
# - npVector
# - vector_trackClones           vector that tracks its clones
# - vector_no_trackClones        vector that does not track its clones
# - vector_trackClones_dual      vector that tracks its clones and has non-default dual
# - vector_no_trackClones_dual   vector that does not track its clones and has non-default dual

class vector_base(Vector):
    def __init__(self, values=None):
        assert isinstance(values, np.ndarray)
        assert values.ndim == 1
        self.values = values
        super().__init__()

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

    def zero(self):
        self.setScalar(0.)

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

    def basis(self, i):
        b = self.clone()
        b.values[:] = 0.
        b.values[i] = 1.
        self._basis = b
        return b


# currently, npVector is implemented like this
class vector_trackClones(vector_base):
    def __init__(self, values=None):
        super().__init__(values)
        self.copies = []

    def __del__(self):
        for copy in self.copies:
            del copy

    def clone(self):
        tmp = type(self)(np.full(self.values.shape, fill_value=np.nan))
        self.copies.append(tmp)
        return tmp


class vector_trackClones_with_dual(vector_trackClones):
    def dual(self):
        return type(self)(self.values.copy())


class vector_no_trackClones(vector_base):
    def __init__(self, values=None):
        super().__init__(values)

    def clone(self):
        tmp = type(self)(np.full(self.values.shape, fill_value=np.nan))
        return tmp


class vector_no_trackClones_with_dual(vector_no_trackClones):
    def dual(self):
        return type(self)(self.values.copy())


class vector_no_trackClones_with_dual_alternative(vector_no_trackClones):
    def dual(self):
        return self.clone()


# We are testing the constructor of Reduced_Objective_SimOpt.
#
# Here is what that does:
#
# stateStore_   = makePtr<VectorController<Real>>();
# adjointStore_ = makePtr<VectorController<Real>>();
# state_        = state->clone(); state_->set(*state);
# adjoint_      = adjoint->clone();
# state_sens_   = state->clone();
# adjoint_sens_ = adjoint->clone();
# dualstate_    = state->dual().clone();
# dualstate1_   = state->dual().clone();
# dualadjoint_  = adjoint->dual().clone();
# dualcontrol_  = control->dual().clone();

# Set up some fake objective and constraint.
# They don't matter for the constructor of Reduced_Objective_SimOpt.
obj = Objective_SimOpt()
constraint = Constraint_SimOpt()


# no segfault
u = npVector(np.ones(3))
z = npVector(np.ones(2))
p = u.dual().clone()
Reduced_Objective_SimOpt(obj, constraint, u, z, p)
print('npVector passed')


# no segfault
u = vector_trackClones(np.ones(3))
z = vector_trackClones(np.ones(2))
p = u.dual().clone()
Reduced_Objective_SimOpt(obj, constraint, u, z, p)
print('vector_trackClones passed')


# no segfault
u = vector_no_trackClones(np.ones(3))
z = vector_no_trackClones(np.ones(2))
p = u.dual().clone()
Reduced_Objective_SimOpt(obj, constraint, u, z, p)
print('vector_no_trackClones passed')


# segfault
u = vector_trackClones_with_dual(np.ones(3))
z = vector_trackClones_with_dual(np.ones(2))
p = u.dual().clone()
Reduced_Objective_SimOpt(obj, constraint, u, z, p)
print('vector_trackClones_with_dual passed')


# segfault
u = vector_no_trackClones_with_dual(np.ones(3))
z = vector_no_trackClones_with_dual(np.ones(2))
p = u.dual().clone()
Reduced_Objective_SimOpt(obj, constraint, u, z, p)
print('vector_no_trackClones_with_dual passed')


# segfault
u = vector_no_trackClones_with_dual_alternative(np.ones(3))
z = vector_no_trackClones_with_dual_alternative(np.ones(2))
p = u.dual().clone()
Reduced_Objective_SimOpt(obj, constraint, u, z, p)
print('vector_no_trackClones_with_dual_alternative passed')

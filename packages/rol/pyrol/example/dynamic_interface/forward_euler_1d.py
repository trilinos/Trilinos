# TO-DO: put in separate file(s) ##############################################

from abc import ABC, abstractmethod


class VectorField(ABC):

    @abstractmethod
    def value(self, u, z, t):
        pass

    @abstractmethod
    def applyJacobian_u(self, jv, v, u, z, t):
        pass

    @abstractmethod
    def applyJacobian_z(self, jv, v, u, z, t):
        pass


class DecayingExponential(VectorField):

    def value(self, u, z, t):
        return z[:] - u[:]

    def applyJacobian_u(self, jv, v, u, z, t):
        jv[:] = -v[:]

    def applyJacobian_z(self, jv, v, u, z, t):
        jv[:] = +v[:]


# TO-DO: put in separate file #################################################

from pyrol.unsupported import DynamicConstraint


class ForwardEuler(DynamicConstraint):

    def __init__(self, f):
        super().__init__()
        self.f = f

    def value(self, c, uo, un, z, ts):
        dt = ts.t[1] - ts.t[0]
        c[:] = un[:] - uo[:] - dt*self.f.value(uo, z, ts.t[0])

    def solve(self, c, uo, un, z, ts):
        dt = ts.t[1] - ts.t[0]
        un[:] = uo[:] + dt*self.f.value(uo, z, ts.t[0])
        self.value(c, uo, un, z, ts)

    def applyJacobian_uo(self, jv, vo, uo, un, z, ts):
        self.f.applyJacobian_u(jv, vo, uo, z, ts.t[0])
        dt = ts.t[1] - ts.t[0]
        jv[:] = - vo[:] - dt*jv[:]

    def applyJacobian_un(self, jv, vn, uo, un, z, ts):
        jv[:] = vn[:]

    def applyJacobian_z(self, jv, vz, uo, un, z, ts):
        self.f.applyJacobian_z(jv, vz, uo, z, ts.t[0])
        dt = ts.t[1] - ts.t[0]
        jv[:] = - dt*jv[:]

    def applyAdjointJacobian_uo(self, ajv, vo, uo, un, z, ts):
        self.applyJacobian_uo(ajv, vo, uo, un, z, ts)

    def applyAdjointJacobian_un(self, ajv, vn, uo, un, z, ts):
        self.applyJacobian_un(ajv, vn, uo, un, z, ts)

    def applyAdjointJacobian_z(self, ajv, vz, uo, un, z, ts):
        self.applyJacobian_z(ajv, vz, uo, un, z, ts)

    def applyInverseJacobian_un(self, ijv, vn, uo, un, z, ts):
        self.applyJacobian_un(ijv, vn, uo, un, z, ts)

    def applyInverseAdjointJacobian_un(self, iajv, vn, uo, un, z, ts):
        self.applyJacobian_un(iajv, vn, uo, un, z, ts)

    # Gauss-Newton for now


# TO-DO: put in separate file #################################################

from pyrol.unsupported import DynamicObjective


class SquaredLoss(DynamicObjective):

    def __init__(self, T, y):
        super().__init__()
        self.T = T
        self.y = y

    def value(self, uo, un, z, ts):
        v = 0
        if ts.t[1] == self.T:
            v = 0.5*np.sum((un[:] - self.y)**2)
        return v

    def gradient_uo(self, g, uo, un, z, ts):
        g[:] = 0

    def gradient_un(self, g, uo, un, z, ts):
        g[:] = 0
        if ts.t[1] == self.T:
            g[:] = un[:] - self.y

    def gradient_z(self, g, uo, un, z, ts):
        g[:] = 0


# TO-DO: put in separate file #################################################

import numpy as np
import matplotlib.pyplot as plt

from pyrol.vectors import NumPyVector
from pyrol import getCout, Problem, Solver

from pyrol.unsupported import TimeStamp
from pyrol.unsupported import ReducedDynamicObjective
from pyrol.unsupported import PartitionedVector
from pyrol.unsupported import SerialConstraint

from pyrol.pyrol.Teuchos import ParameterList
import pyrol


def as_numpy(v):
    # v : ROL:PartitionedVector
    n = v.dimension()
    w = np.full(n, np.nan)
    for k in range(n):
        w[k] = v[k][:][0]
    return w


def get_trajectory(constraint, u0, z, timestamps):
    serial_constraint = SerialConstraint(constraint, u0, timestamps)
    n = z.dimension()
    c = PartitionedVector.create(NumPyVector(np.array([np.nan])), n)
    u = PartitionedVector.create(NumPyVector(np.array([np.nan])), n)
    serial_constraint.solve(c, u, z, np.nan)
    t = np.full(n, np.nan)
    for k in range(n):
        t[k] = timestamps[k].t[0]
    u = np.hstack([u0[:][0], as_numpy(u)[:-1]])
    return t, u


def plot_trajectory(constraint, u0, z, timestamps):
    t, u = get_trajectory(constraint, u0, z, timestamps)
    plt.plot(t, np.exp(-t), label='continuous')
    plt.plot(t, u, '.-', label='numerical solution')
    plt.grid()
    plt.title(f'# of timesteps = {z.dimension()}')
    plt.ylabel('state')
    plt.xlabel('t')
    plt.legend()
    plt.show()


def optimize(n):

    T  = 1   # end time

    # configure timestamps
    dt = T/n
    timestamps = pyrol.pyrol.std.vector_ROL_TimeStamp_double_t()
    for k in range(n):
        ts = TimeStamp()
        ts.t.push_back(k*dt)
        ts.t.push_back((k + 1)*dt)
        timestamps.push_back(ts)

    # configure DynamicConstraint
    f = DecayingExponential()
    y = np.exp(-1)  # the output we want to match (an everywhere zero control)
    constraint = ForwardEuler(f)

    # configure DynamicObjective
    objective = SquaredLoss(T, y)

    parameters = ParameterList()
    cout = getCout()

    # configure ReducedDynamicObjective
    u0 = NumPyVector(np.ones(1))
    zk = NumPyVector(np.zeros(1))
    ck = NumPyVector(np.zeros(1))
    reduced_objective = ReducedDynamicObjective(objective, constraint, u0, zk, ck, timestamps, parameters, cout)

    # Problem
    z = PartitionedVector.create(zk, n)
    problem = Problem(reduced_objective, z)
    # problem.check(True, cout)

    parameters['General'] =  ParameterList()
    parameters['General']['Output Level'] = 1
    parameters['Step'] = ParameterList()
    parameters['Step']['Trust Region'] = ParameterList()
    parameters['Step']['Trust Region']['Subproblem Solver'] = 'Dogleg'

    # Solver
    solver = Solver(problem, parameters)
    solver.solve(cout)

    return constraint, u0, z, timestamps


def main():
    n = 100  # number of time steps
    constraint, u0, z, timestamps = optimize(n)
    print(f'\n||z||_inf = {np.linalg.norm(as_numpy(z), np.inf)}\n')
    plot_trajectory(constraint, u0, z, timestamps)


if __name__ == '__main__':
    main()

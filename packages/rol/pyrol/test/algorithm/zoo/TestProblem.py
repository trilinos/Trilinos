# compare with src/zoo/testproblems/ROL_TestProblem.hpp

import numpy as np

from abc import ABC, abstractmethod
from pyrol import Problem, Bounds
from pyrol.vectors import NumPyVector


class TestProblem(ABC):

    @abstractmethod
    def getObjective(self):
        pass

    @abstractmethod 
    def getInitialGuess(self):
        pass

    @abstractmethod
    def getSolution(self, i = 0):
        pass

    def getNumSolutions(self):
       return 1

    def getBoundConstraint(self):
       return None

    def getLinearEqualityConstraint(self):
       return None

    def getLinearEqualityMultiplier(self):
       return None
    
    def getEqualityConstraint(self):
       return None

    def getEqualityMultiplier(self):
       return None

    def getLinearInequalityConstraint(self):
       return None

    def getLinearInequalityMultiplier(self):
       return None

    def getInequalityConstraint(self):
       return None

    def getInequalityMultiplier(self):
       return None

    def get(self):
        x0 = self.getInitialGuess()
        x = []
        for i in range(self.getNumSolutions()):
            x.append(self.getSolution(i))
        problem = Problem(self.getObjective(), x0)

        bcon = self.getBoundConstraint()
        if bcon is not None:
            problem.addBoundConstraint(bcon)

        econ = self.getEqualityConstraint()
        if econ is not None:
            emul = self.getEqualityMultiplier()
            assert emul is not None
            problem.addConstraint("equality constraint", econ, emul)

        linear_econ = self.getLinearEqualityConstraint()
        if linear_econ is not None:
            linear_emul = self.getLinearEqualityMultiplier()
            assert linear_emul is not None
            problem.addLinearConstraint("linear equality constraint", 
                                        linear_econ, linear_emul)

        icon = self.getInequalityConstraint()
        if icon is not None: 
            imul = self.getInequalityMultiplier()
            assert imul is not None
            b = NumPyVector(np.zeros(imul.array.size))
            bcon = Bounds(b, True)
            problem.addConstraint("inequality constaint", icon, imul, bcon)

        linear_icon = self.getLinearInequalityConstraint()
        if linear_icon is not None:
            linear_imul = self.getLinearInequalityMultiplier()
            assert linear_imul is not None
            b = NumPyVector(np.zeros(linear_imul.array.size))
            bcon = Bounds(b, True)
            problem.addLinearConstraint("linear inquality constraint", 
                                        linear_icon, linear_imul, bcon)

        return problem, x0, x
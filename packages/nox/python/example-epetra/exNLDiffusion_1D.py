#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#                 PyTrilinos.NOX: Python Interface to NOX
#                   Copyright (2005) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
# @HEADER

# Imports.  Users importing an installed version of PyTrilinos should use the
# "from PyTrilinos import ..." syntax.  Here, the setpath module adds the build
# directory, including "PyTrilinos", to the front of the search path.  We thus
# use "import ..." for Trilinos modules.  This prevents us from accidentally
# picking up a system-installed version and ensures that we are testing the
# build module.
import sys

try:
    import setpath
    import Epetra
    import EpetraExt
    import NOX
except ImportError:
    from PyTrilinos import Epetra, EpetraExt, NOX
    print >>sys.stderr, "Using system-installed Epetra, EpetraExt, NOX"

from Numeric import *

########################################################################

class Problem:
    def __init__(self, k, size, comm):
        self.__k    = k
        self.__size = size
        self.__comm = comm
        self.__mesh = 1.0 / (size-1.0)
        self.__grid = arange(size) * self.__mesh
        self.__map  = Epetra.Map(size,0,comm)
        self.__soln = Epetra.Vector(self.__map)
        self.__x    = zeros(size,Float64)
        self.__rhs  = zeros(size,Float64)
        self.computeGraph()
    def computeGraph(self):
        self.__graph = Epetra.CrsGraph(Epetra.Copy, self.__map, 3)
        for lrid in range(self.__graph.NumMyRows()):
            grid = self.__graph.GRID(lrid)
            if   grid == 0            : indices = [0,1]
            elif grid == self.__size-1: indices = [grid-1,grid]
            else                      : indices = [grid-1,grid,grid+1]
            self.__graph.InsertGlobalIndices(grid,len(indices),indices)
        self.__graph.FillComplete()
    def setInterface(self, interface):
        self.__interface = interface
    def getSolution(self):
        return self.__soln
    def getGraph(self):
        return self.__graph
    def computeF(self):
        rhs = self.__rhs
        x   = self.__x
        self.__interface.unloadX(x)
        # BC: x(0.0) = 1.0
        rhs[0] = x[0] - 1.0
        # Interior:  x" + k x^2 = 0
        for i in range(1,self.__size-1):
            rhs[i] = self.__k * x[i] * x[i] + \
                     (x[i-1] - 2*x[i] + x[i+1]) / (self.__mesh*self.__mesh)
        # BC: x(1.0) = 1.0
        rhs[self.__size-1] = x[self.__size-1] - 1.0
        self.__interface.loadRHS(rhs)
        return True
    def computeJacobian(self, x, m):
        return False
    def computePrecMatrix(self, x, m):
        return False
    def computePreconditioner(self, x, op):
        return False

########################################################################

def main():
    # Create a communicator
    comm    = Epetra.PyComm()
    myPID   = comm.MyPID()
    numProc = comm.NumProc()

    # Get the problem size from the command line
    if (len(sys.argv) > 2):
        print "usage:", sys.argv[0], "problemSize"
        sys.exit(1)
    if (len(sys.argv) == 1):
        probSize = 11 * comm.NumProc()
    else:
        probSize = int(sys.argv[1])
    if myPID == 0: print "probSize =", probSize

    # Create an instance of the nonlinear Problem class
    # and initialize the first guess
    problem = Problem(1.0, probSize, comm)
    soln = problem.getSolution()
    soln.PutScalar(1.0)

    # Create the top level parameter list
    nlParams     = NOX.Parameter.List()

    # Set the nonlinear solver method
    nlParams.setParameter("Nonlinear Solver", "Line Search Based")

    # Set the printing parameters in the "Printing" sublist
    printParams  = nlParams.sublist("Printing")
    printParams.setParameter("MyPID",              myPID) 
    printParams.setParameter("Output Precision",   3    )
    printParams.setParameter("Output Processor",   0    )
    printParams.setParameter("Output Information", 0    )

    # Sublist for line search
    searchParams = nlParams.sublist("Line Search")
    searchParams.setParameter("Method", "Full Step")

    # Sublist for direction
    dirParams    = nlParams.sublist("Direction")
    dirParams.setParameter("Method", "Newton")

    # Sublist for Newton parameters
    newtonParams = dirParams.sublist("Newton")
    newtonParams.setParameter("Forcing Term Method", "Constant")

    # Sublist for linear solver
    lsParams     = newtonParams.sublist("Linear Solver")
    lsParams.setParameter("Aztec Solver",     "GMRES"                   )
    lsParams.setParameter("Max Iterations",   800                       )
    lsParams.setParameter("Tolerance",        1e-4                      )
    lsParams.setParameter("Scaling",          "None"                    )
    lsParams.setParameter("Preconditioning",  "AztecOO: Jacobian Matrix")

    # Create the interface between the Problem and the NOX nonlinear solver
    interface = NOX.Epetra.PyInterface(problem)
    problem.setInterface(interface)

    # Generate a color map from the problem graph
    mapColoring   = EpetraExt.CrsGraph_MapColoring(False)
    colorMap      = mapColoring(problem.getGraph())
#    print "colorMap = {", colorMap, "\n}"
    colorMapIndex = EpetraExt.CrsGraph_MapColoringIndex(colorMap)
    columns       = colorMapIndex(problem.getGraph())

    # Create the finite difference coloring object
    fdc = NOX.Epetra.FiniteDifferenceColoring(interface, soln,    \
                                              problem.getGraph(), \
                                              colorMap, columns)

    # Create the group
    group = NOX.Epetra.Group(printParams, lsParams, interface, soln, fdc)

    # Create the convergence tests
    absresid  = NOX.StatusTest.NormF(1.0e-6)
    relresid  = NOX.StatusTest.NormF(group, 1.0e-6)
    update    = NOX.StatusTest.NormUpdate(1.0e-5)
    wrms      = NOX.StatusTest.NormWRMS(1.0e-2, 1.0e-8)
    maxiters  = NOX.StatusTest.MaxIters(20)
    converged = NOX.StatusTest.Combo(NOX.StatusTest.Combo.AND)
    if myPID == 0: print "Adding absresid"
    converged.addStatusTest(absresid)
    if myPID == 0: print "Adding relresid"
    converged.addStatusTest(relresid)
    if myPID == 0: print "Adding wrms"
    converged.addStatusTest(wrms    )
    if myPID == 0: print "Adding update"
    converged.addStatusTest(update  )
    if myPID == 0: print "Creating combo"
    combo     = NOX.StatusTest.Combo(NOX.StatusTest.Combo.OR)
    if myPID == 0: print "Adding converged"
    combo.addStatusTest(converged)
    if myPID == 0: print "Adding maxiters"
    combo.addStatusTest(maxiters)
    if myPID == 0: print "Done building combo"

    # Create the method
    solver = NOX.Solver.Manager(group, combo, nlParams)

    # Solve the problem
    status = solver.solve()

    # Print the status test result
    if myPID == 0:
        print "\nStatus test:"
        if (status == NOX.StatusTest.Converged):
            print "Nonlinear solver converged!"
        elif (status == NOX.StatusTest.Failed):
            print "Nonlinear solver failed to converge!"
        elif (status == NOX.StatusTest.Unconverged):
            print "Nonlinear solver neither failed nor converged!"
        elif (status == NOX.StatusTest.Unevaluated):
            print "Nonlinear solver unevaluated."

    # Print the final output parameters
    nlOutputParams = nlParams.sublist("Output")
    lsOutputParams = lsParams.sublist("Output")
    nlIts = nlOutputParams.getParameter("Nonlinear Iterations")
    lsIts = lsOutputParams.getParameter("Total Number of Linear Iterations")
    resid = nlOutputParams.getParameter("2-Norm of Residual")
    if myPID == 0:
        print "\n", \
              "Final Parameters\n", \
              "****************\n", \
              "Nonlinear iterations:   ", nlIts, "\n", \
              "Total linear iterations:", lsIts, "\n", \
              "Final residual:         ", resid

if (__name__ == "__main__"):
    main()

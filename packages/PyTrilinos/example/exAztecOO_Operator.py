#! /usr/bin/env python

# @HEADER
# ************************************************************************
# 
#                PyTrilinos: Python Interface to Trilinos
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
# Questions? Contact Bill Spotz (wfspotz@sandia.gov) 
# 
# ************************************************************************
# @HEADER

################################################################################
#
# Example of how to define an Epetra.Operator derived class in python
#
# An Epetra_Operator is an object that can be applied to vectors to return a
# vector, but has no structure -- that is, there is no "row" or "graph".
# Basically, an Epetra_Operator is defined by two methods: Apply() and
# ApplyInverse().  In this case, the operator represents a 1D Laplacian on a
# structure grid.  The linear system is solved using AztecOO and no
# preconditioning (since algebraic preconditioners are hard to define for pure
# Operator classes).
#
# Based on a script originally written by Marzio Sala.  Updated by Bill Spotz.
#
################################################################################

try:
    import setpath
    import Epetra
    import AztecOO
except:
    from PyTrilinos import Epetra, AztecOO
    print "Using installed versions of Epetra, AztecOO"

################################################################################

class Laplace1D_Operator(Epetra.Operator):

    def __init__(self, n, comm=None):
        """
        __init__(self, n) -> Laplace1D_Operator (with an Epetra.PyComm() communicator)
        __init__(self, n, comm) -> Laplace1D_Operator (with given communicator)
        """
        # Initialize the base class.  This is REQUIRED
        Epetra.Operator.__init__(self)
        # Determine the communicator
        if comm is None:
            self.__comm = Epetra.PyComm()
        else:
            self.__comm = comm
        # Default indexes
        self.__y0 =  1
        self.__y1 = -1
        # Create the range map
        self.__rangeMap = Epetra.Map(n,0,self.__comm)
        # Create the domain map
        myIndexes = list(self.__rangeMap.MyGlobalElements())
        if self.__comm.MyPID() > 0:
            self.__y0 = None    # Equivalent to index 0
            myIndexes.insert(0,myIndexes[0]-1)
        if self.__comm.MyPID() < self.__comm.NumProc()-1:
            self.__y1 = None    # Equivalent to last index
            myIndexes.append(myIndexes[-1]+1)
        self.__domainMap = Epetra.Map(-1, myIndexes, 0, self.__comm)
        # Store a label for the operator
        self.__label = "1D Laplace Operator"
        # Transpose flag
        self.__useTranspose = False

    def __str__(self):
        "Return the operator's label"
        return self.__label

    def Label(self):
        "Required implementation of Epetra.Operator class"
        return self.__label

    def OperatorDomainMap(self):
        "Required implementation of Epetra.Operator class"
        return self.__domainMap

    def OperatorRangeMap(self):
        "Required implementation of Epetra.Operator class"
        return self.__rangeMap

    def Comm(self):
        "Required implementation of Epetra.Operator class"
        return self.__comm

    def Apply(self,x,y):
        """
        Required implementation of Epetra.Operator class.  This method will be
        called by the AztecOO solver in order to compute y = Ax, where A is this
        operator.
        """
        try:
            # Apply operator to interior points
            y[:,self.__y0:self.__y1] = 2.0 * x[:,1:-1] - x[:,:-2] - x[:,2:]
            # Apply left boundary condition
            if self.__comm.MyPID() == 0:
                y[:,0] = x[:,0]
            # Apply right boundary condition
            if self.__comm.MyPID() == self.__comm.NumProc() - 1:
                y[:,-1] = x[:,-1]
        except Exception, e:
            print e
            return -1

        return 0

    def ApplyInverse(self):
        "Required implementation of Epetra.Operator class"
        return -1

    def HasNormInf(self):
        "Required implementation of Epetra.Operator class"
        return True

    def NormInf(self):
        "Required implementation of Epetra.Operator class"
        return 4.0

    def SetUseTranspose(self, useTranspose):
        "Required implementation of Epetra.Operator class"
        self.__useTranspose = bool(useTranspose)

    def UseTranspose(self):
        "Required implementation of Epetra.Operator class"
        return self.__useTranspose

################################################################################

def main():

    # Problem initialization
    n     = 100
    bc0   = 0.0
    bc1   = 1.0
    tol   = 1.0e-5
    comm  = Epetra.PyComm()
    lap1d = Laplace1D_Operator(n, comm)

    # Create solution and RHS vectors
    x = Epetra.Vector(lap1d.OperatorDomainMap())
    b = Epetra.Vector(lap1d.OperatorRangeMap() )

    # Initialize vectors: x will be a straight line between its boundary values,
    # and b=1, with its boundary values equal to x on the boundaries
    x[:] = bc0 + (bc1-bc0) * (x.Map().MyGlobalElements() / (n-1.0))
    b.PutScalar(1.0)
    if comm.MyPID() == 0:
        b[0] = bc0
    if comm.MyPID() == comm.NumProc()-1:
        b[-1] = bc1

    # Build the linear system solver
    problem = Epetra.LinearProblem(lap1d, x, b)
    solver  = AztecOO.AztecOO(problem)
    solver.SetParameters({"Solver"  : "CG",
                          "Precond" : "None",
                          "Output"  : 16    })

    # Solve the problem
    solver.Iterate(n, tol)
    if comm.MyPID() == 0:
        if solver.ScaledResidual() < tol: print "End Result: TEST PASSED"
        else:                             print "End Result: TEST FAILED"

################################################################################

if __name__ == "__main__":
    main()

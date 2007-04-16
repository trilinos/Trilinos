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

######################################################################
#
# Example of how to define an Epetra.RowMatrix derived class in python
# 
# This example shows how to derive from the class Epetra.RowMatrix in python.
# The main procedure is as follows: - create a python class derived from
# Epetra.RowMatrix -- in this case a 1D Laplacian matrix.  Define all the
# methods as done in this file.  The most important methods are probably
# Apply(), Multiply(), and ExtractMyRowCopy() (note that the return signature
# for NumMyRowEntries() and ExtractMyRowCopy() are different for python than for
# C++).  Some methods do not need to be implemented; in this case they simply
# return -1.  You may now create an instance of your derived class and supply it
# to any Trilinos solver that can use it (in this case AztecOO).
#
# Based on a script originally written by Marzio Sala.  Updated by Bill Spotz.
#
######################################################################

try:
    import setpath
    import Epetra
    import AztecOO
except:
    from PyTrilinos import Epetra, AztecOO
    print "Using installed version of Epetra, AztecOO"

import numpy

################################################################################

class Laplace1D_Matrix(Epetra.RowMatrix):

    def __init__(self, n, comm=None):
        """
        __init__(self, n) -> Laplace1D_Matrix (with and Epetra.PyComm() communicator)
        __init__(self, n, comm) -> Laplace1D_Matrix (with given communicator)
        """
        # Initialize the base class.  This is REQUIRED
        Epetra.RowMatrix.__init__(self)
        # Determine the communicator
        if comm is None:
            self.__comm = Epetra.PyComm()
        else:
            self.__comm = comm
        # Default indexes
        self.__y0 =  1
        self.__y1 = -1
        # Create the row map
        self.__rowMap = Epetra.Map(n, 0, self.__comm)
        # Create the col map
        myIndexes = list(self.__rowMap.MyGlobalElements())
        if self.__comm.MyPID() > 0:
            self.__y0 = None    # Equivalent to index 0
            myIndexes.insert(0,myIndexes[0]-1)
        if self.__comm.MyPID() < self.__comm.NumProc()-1:
            self.__y1 = None    # Equivalent to last index
            myIndexes.append(myIndexes[-1]+1)
        self.__colMap = Epetra.Map(-1, myIndexes, 0, self.__comm)
        # Store a label for the row matrix
        self.__label = "1D Laplace Row Matrix"
        # Store the matrix properties
        self.__numRows = n
        self.__numCols = n
        self.__useTranspose = False

    def __str__(self):
        "Return the row matrix label"
        return self.__label

    def Map(self):
        "Required implementation of Epetra.SrcDistObject class"
        return self.__rowMap

    def SetUseTranspose(self, useTranspose):
        "Required implementation of Epetra.Operator class"
        self.__useTranspose = bool(useTranspose)

    def UseTranspose(self):
        "Required implementation of Epetra.Operator class"
        return self.__useTranspose

    def Apply(self, LHS, RHS):
        "Required implementation of Epetra.Operator class"
        return self.Multiply(self.__useTranspose, LHS, RHS)

    def ApplyInverse(self):
        "Required implementation of Epetra.Operator class"
        return -2

    def HasNormInf(self):
        "Required implementation of Epetra.Operator class"
        return True

    def NormInf(self):
        "Required implementation of Epetra.Operator class"
        return 4.0

    def Label(self):
        "Required implementation of Epetra.Operator class"
        return self.__label

    def Comm(self):
        "Required implementation of Epetra.Operator class"
        return self.__comm

    def OperatorDomainMap(self):
        "Required implementation of Epetra.Operator class"
        return self.__colMap

    def OperatorRangeMap(self):
        "Required implementation of Epetra.Operator class"
        return self.__rowMap

    def NumMyRowEntries(self, MyRow, NumEntries):
        """
        Required implementation of Epetra.RowMatrix class.  In C++, NumEntries
        is an int& argument intended as output.  When called via callbacks from
        C++, this int& is converted to a numpy array of length 1 so that it can
        be altered in-place via NumEntries[0] = ...
        """
        globalRow = self.__rowMap.GID(MyRow)
        if globalRow == 0 or globalRow == self.__numRows-1:
            NumEntries[0] = 1
        else:
            NumEntries[0] = 3
        return 0

    def MaxNumEntries(self):
        "Required implementation of Epetra.RowMatrix class"
        return 3

    def ExtractMyRowCopy(self, MyRow, Length, NumEntries, Values, Indices):
        """
        Required implementation of Epetra.RowMatrix class.  In C++, NumEntries,
        Values, and Indices are all output arguments.  When called via callbacks
        from C++, these arguments are converted to numpy arrays so that we can
        manipulate the data in-place.  NumEntries is a scalar in C++, but must
        be accessed as NumEntries[0] in python.
        """
        globalRow = self.__rowMap.GID(MyRow)
        if globalRow == 0 or globalRow == self.__numRows-1:
            if (Length < 1):
                return -1
            NumEntries[0] = 1
            Values[0]     = 1.0
            Indices[0]    = MyRow
        else:
            if (Length < 3):
                return -1
            NumEntries[0] = 3
            Values[:3]    = [   -1.0,   2.0,    -1.0]
            Indices[:3]   = [MyRow-1, MyRow, MyRow+1]
        return 0

    def ExtractDiagonalCopy(self, Diagonal):
        "Required implementation of Epetra.RowMatrix class"
        Diagonal.PutScalar(2.0)
        myPID = self.__comm.MyPID()
        if myPID == 0:                       Diagonal[ 0] = 1.0
        if myPID == self.__comm.NumProc()-1: Diagonal[-1] = 1.0
        return 0

    def Multiply(self, UseTranspose, x, y):
        "Required implementation of Epetra.RowMatrix class"
        try:
            # Apply operator to interior points
            y[:,self.__y0:self.__y1] = 2.0 * x[:,1:-1] - x[:,:-2] - x[:,2:]
            # Apply left boundary condition
            if self.__comm.MyPID() == 0:
                y[:,0] = x[:,0]
            # Apply right boundary condition
            if self.__comm.MyPID() == self.__comm.NumProc() - 1:
                y[:,-1] = x[:,-1]
            return 0
        except Exception, e:
            print e
            return -1

    def Solve(self, upper, trans, unitDiagonal, x, y):
        "Required implementation of Epetra.RowMatrix class"
        return -1

    def InvRowSums(self, x):
        "Required implementation of Epetra.RowMatrix class"
        return -1

    def LeftScale(self, x):
        "Required implementation of Epetra.RowMatrix class"
        return -1

    def InvColSums(self, x):
        "Required implementation of Epetra.RowMatrix class"
        return -1

    def RightScale(self, x):
        "Required implementation of Epetra.RowMatrix class"
        return -1

    def Filled(self):
        "Required implementation of Epetra.RowMatrix class"
        return True

    def NormOne(self):
        "Required implementation of Epetra.RowMatrix class"
        return 4.0

    def NumGlobalNonzeros(self):
        "Required implementation of Epetra.RowMatrix class"
        return 3 * self.__numRows - 2

    def NumGlobalRows(self):
        "Required implementation of Epetra.RowMatrix class"
        return self.__numRows

    def NumGlobalCols(self):
        "Required implementation of Epetra.RowMatrix class"
        return self.__numCols

    def NumGlobalDiagonals(self):
        "Required implementation of Epetra.RowMatrix class"
        return self.__numRows

    def NumMyNonzeros(self):
        "Required implementation of Epetra.RowMatrix class"
        return 3 * self.__numRows - 2

    def NumMyRows(self):
        "Required implementation of Epetra.RowMatrix class"
        return self.__numRows

    def NumMyCols(self):
        "Required implementation of Epetra.RowMatrix class"
        return self.__numCols

    def NumMyDiagonals(self):
        "Required implementation of Epetra.RowMatrix class"
        return self.__numRows

    def LowerTriangular(self):
        "Required implementation of Epetra.RowMatrix class"
        return False

    def UpperTriangular(self):
        "Required implementation of Epetra.RowMatrix class"
        return False

    def RowMatrixRowMap(self):
        "Required implementation of Epetra.RowMatrix class"
        return self.__rowMap

    def RowMatrixColMap(self):
        "Required implementation of Epetra.RowMatrix class"
        return self.__colMap

    def RowMatrixImporter(self):
        "Required implementation of Epetra.RowMatrix class"
        return Epetra.Import(self.__rowMap, self.__colMap)

################################################################################

def main():

    # Problem initialization
    n     = 100
    bc0   = 0.0
    bc1   = 1.0
    tol   = 1.0e-5
    comm  = Epetra.PyComm()
    lap1d = Laplace1D_Matrix(n, comm)

    # Create solution and RHS vectors
    x = Epetra.Vector(lap1d.RowMatrixColMap())
    b = Epetra.Vector(lap1d.RowMatrixRowMap())

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
                          "Precond" : "Jacobi",
                          "Output"  : 16      })

    # Solve the problem
    solver.Iterate(5*n, tol)
    if comm.MyPID() == 0:
        if solver.ScaledResidual() < tol: print "End Result: TEST PASSED"
        else:                             print "End Result: TEST FAILED"

################################################################################

if __name__ == "__main__":
    main()

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
# Example of 1D Laplacian on one processor.
# 
# This example shows how to derive class Epetra_RowMatrix in Python.
# For the sake of simplicity, the example runs only with one processor.
# The main procedure is as follows:
# - create a Python class derived from Epetra.RowMatrix
# - define all the methods as done in this file. Probably, the most
#   important methods are Apply, Multiply, and ExtractMyRowCopy. Some
#   methods do not need to be implemented; in this case simply returns
#   `False.'
# - declare the object and use it in all Trilions modules that accept
#   Epetra_RowMatrix's.
#
# Note: You might need to be able to allocate/use integer and double
#       pointers or references. Some tools are available for that:
#       use for example the DArray and IArray classes of the Epetra module.
#
# \author Marzio Sala, ETHZ/COLAB
#
# \date Last updated on 04-Dec-05
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

class Laplace1D(Epetra.RowMatrix):

    def __init__(self, n, comm):
        Epetra.RowMatrix.__init__(self)
        self.__useTranspose = False
        self.__comm         = comm
        self.__numRows      = n
        self.__numCols      = n
        self.__rowMap       = Epetra.Map(self.__numRows, 0, comm)
        self.__colMap       = Epetra.Map(self.__numCols, 0, comm)
        self.__numEntries   = Epetra.IArray(1)
        self.__indices      = Epetra.IArray(3)
        self.__values       = Epetra.DArray(3)
        #self.__numEntries   = numpy.array([0])
        #self.__indices      = numpy.array([0,0,0])
        #self.__values       = numpy.array([0.0,0.0,0.0])

    def __str__(self):
        return "Laplace1D"

    ###############################################################################
    # The following method is an implementation of the Epetra.SrcDistObject class #
    ###############################################################################

    def Map(self):
        return self.__rowMap

    ##########################################################################
    # The following methods are implementations of the Epetra.Operator class #
    ##########################################################################

    def SetUseTranspose(self, useTranspose):
        self.__useTranspose = bool(useTranspose)

    def UseTranspose(self):
        return self.__useTranspose

    def Apply(self, LHS, RHS):
        return self.Multiply(self.__useTranspose, LHS, RHS)

    def ApplyInverse(self):
        return -2

    def HasNormInf(self):
        return True

    def NormInf(self):
        return 4.0

    def Label(self):
        return "Laplace1D"

    def Comm(self):
        return self.__comm

    def OperatorDomainMap(self):
        return self.__rowMap

    def OperatorRangeMap(self):
        return self.__rowMap

    ###########################################################################
    # The following methods are implementations of the Epetra.RowMatrix class #
    ###########################################################################

    def NumMyRowEntries(self, MyRow, n):
        NumEntries = Epetra.IArray(n, 1)
        if MyRow == 0 or MyRow == self.__numRows:
            NumEntries[0] = 2
        else:
            NumEntries[0] = 3
        return 0

    def MaxNumEntries(self):
        return 3

    # Input to this function one has an int* pointer (NumEntries),
    # a double* pointer (Values) and an int* pointer(Indices); these
    # pointers are wrapped as IArray's and DArray's
    def ExtractMyRowCopy(self, MyRow, Length, n, vals, ind):
        if (Length < 3):
            return -1
        NumEntries = Epetra.IArray(n,    1     )
        Values     = Epetra.DArray(vals, Length)
        Indices    = Epetra.IArray(ind,  Length)
        n = self.__numRows
        if MyRow == 0:
            Indices[0] = 0
            Indices[1] = 1
            Values[0] = 2.0
            Values[1] = -1.0
            NumEntries[0] = 2
        elif MyRow == n - 1:
            Indices[0] = n - 1; Indices[1] = n - 2
            Values[0] = 2.0; Values[1] = -1.0
            NumEntries[0] = 2
        else:
            Indices[0] = MyRow; Indices[1] = MyRow - 1; Indices[2] = MyRow + 1
            Values[0] = 2.0; Values[1] = -1.0; Values[2] = -1.0
            NumEntries[0] = 3
        return 0

    def ExtractDiagonalCopy(self, Diagonal):
        Diagonal.PutScalar(2.0)
        return 0

    def Multiply(self, UseTranspose, LHS, RHS):
        Indices      = self.__indices
        Values       = self.__values
        NumEntries   = self.__numEntries

        n = RHS.MyLength()
        if LHS.NumVectors() != 1:
            print "this Apply() function has been implemented for a single vector"
            return -1

        # I need to wrap the Values() array as DArray; if I use the bracket
        # operator the code crashes...
        LHS_V = Epetra.DArray(LHS.Values(), n)
        RHS_V = Epetra.DArray(RHS.Values(), n)

        for i in xrange(self.__numRows):
            ierr = self.ExtractMyRowCopy(i, 5, NumEntries.Values(),
                                         Values.Values(), Indices.Values())
            total = 0.0
            for j in xrange(NumEntries[0]):
                total = total + LHS_V[Indices[j]] * Values[j]
            RHS_V[i] = total  

        return 0

    def Solve(self, upper, trans, unitDiagonal, x, y):
        return -1

    def InvRowSums(self, x):
        return -1

    def LeftScale(self, x):
        return -1

    def InvColSums(self, x):
        return -1

    def RightScale(self, x):
        return -1

    def Filled(self):
        return True

    def NormOne(self):
        return 4.0

    def NumGlobalNonzeros(self):
        return 3 * self.__numRows - 2

    def NumGlobalRows(self):
        return self.__numRows

    def NumGlobalCols(self):
        return self.__numCols

    def NumGlobalDiagonals(self):
        return self.__numRows

    def NumMyNonzeros(self):
        return 3 * self.__numRows - 2

    def NumMyRows(self):
        return self.__numRows

    def NumMyCols(self):
        return self.__numCols

    def NumMyDiagonals(self):
        return self.__numRows

    def LowerTriangular(self):
        return False

    def UpperTriangular(self):
        return False

    def RowMatrixRowMap(self):
        return self.__rowMap

    def RowMatrixColMap(self):
        return self.__colMap

    def RowMatrixImporter(self):
        return Epetra.Import(self.__rowMap, self.__rowMap)

################################################################################

def main():

    # Problem initialization
    n = 10
    comm = Epetra.PyComm()
    if comm.NumProc() != 1:
        print "This example is only serial, sorry"
        return

    # Create a 1D Laplacian linear system
    matrix = Laplace1D(n, comm)
    lhs    = Epetra.Vector(matrix.Map())
    rhs    = Epetra.Vector(matrix.Map())
    lhs.PutScalar(0.0)
    rhs.PutScalar(1.0)
    solver = AztecOO.AztecOO(matrix, lhs, rhs)

    # Set the AztecOO options
    solver.SetAztecOption(AztecOO.AZ_solver,AztecOO.AZ_gmres)
    solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_Jacobi)
    solver.SetAztecOption(AztecOO.AZ_output, 16)

    # Solve the problem
    solver.Iterate(1550, 1e-5)

    if comm.MyPID() == 0: print "End Result: TEST PASSED"

################################################################################

if __name__ == "__main__":
    main()

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
# Example of how to define an Epetra_Operator derived class in Python
#
# An Epetra_Operator is an object that can be applied to vectors to
# return vector, but has no structure -- that is, there is no "row"
# or "graph". Basically, an Epetra_Operator is defined by two methods:
# Apply() and ApplyInverse(). In this case, the operator represents 
# a 1D Laplacian on a structure grid. The linear system is solved using
# AztecOO and no preconditioning. This example is only serial, but it 
# can be easily extended to parallel.
#
# \authors Marzio Sala, 9214
#          Bill Spotz,  1433
#
# \date Last updated on 22 Feb 2006
#
######################################################################

try:
    import setpath
    import Epetra
    import AztecOO
except:
    from PyTrilinos import Epetra, AztecOO
    print "Using installed versions of Epetra, AztecOO"

################################################################################

class MyOperator(Epetra.Operator):

    def __init__(self, map):
        Epetra.Operator.__init__(self)
        self.__comm            = map.Comm()
        self.__map             = map
        self.__label           = "1D Laplace"
        self.__setUseTranspose = False

    def __str__(self):
        return self.__label

    def Label(self):
        return self.__label

    def Apply(self,x,y):
        # Arguments x and y will be Epetra.Epetra_MultiVectors, i.e. raw wrappers to
        # the C++ class.  In order to have nice python indexing, we need to convert
        # them to Epetra.MultiVectors (or, in this case, Epetra.Vectors).
        lhs = Epetra.Vector(Epetra.View,x,0)
        rhs = Epetra.Vector(Epetra.View,y,0)
        rhs[0]    = 2.0 * lhs[0]    - lhs[1]
        rhs[1:-1] = 2.0 * lhs[1:-1] - lhs[:-2] - lhs[2:]
        rhs[-1]   = 2.0 * lhs[-1]   - lhs[-2]

        return 0

    def ApplyInverse(self):
        return -1

    def OperatorDomainMap(self):
        return self.__map

    def OperatorRangeMap(self):
        return self.__map

    def Map(self):
        return self.__map

    def Comm(self):
        return self.__comm

    def HasNormInf(self):
        return True

    def NormInf(self):
        return 4.0

    def SetUseTranspose(self, useTranspose):
        self.__useTranspose = bool(useTranspose)

    def UseTranspose(self):
        return self.__useTranspose

################################################################################

def main():

    n    = 100
    comm = Epetra.PyComm()
    if comm.NumProc() != 1:
        print "This example is only serial, sorry"
        return
    map = Epetra.Map(n, 0, comm)
    op  = MyOperator(map)

    print op.Label()
    print op.__class__.__bases__
    lhs = Epetra.Vector(map)
    rhs = Epetra.Vector(map)

    rhs.PutScalar(1.0)
    lhs.PutScalar(0.0)

    Problem = Epetra.LinearProblem(op, lhs, rhs)
    print "Creating Solver"
    Solver  = AztecOO.AztecOO(Problem)
    print "Solver Created"

    Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg)
    Solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_none)
    Solver.SetAztecOption(AztecOO.AZ_output, 16)
    Solver.Iterate(1550, 1e-5)

    if comm.MyPID() == 0: print "End Result: TEST PASSED"

################################################################################

if __name__ == "__main__":
    main()

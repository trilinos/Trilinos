#! /usr/bin/env python

# @header
# ************************************************************************
#
#              PyTrilinos.Amesos: Python interface to Amesos
#                   Copyright (2005) Sandia Corporation
#
# Under terms of contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# license, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# without any warranty; without even the implied warranty of
# merchantability or fitness for a particular purpose.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? contact Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
# @header

"""
Usage is: ./exSolvers.py <solver-type>
    where <solver-type> can be:
        - Amesos_Lapack (DEFAULT)
        - Amesos_Klu
        - Amesos_Umfpack
        - Amesos_Superlu
        - Amesos_Superludist
        - Amesos_Dscpack
        - Amesos_Mumps
"""

# System import
import sys

# PyTrilinos imports
#try:
import setpath
#except ImportError:
#  from PyTrilinos import Epetra, Amesos
#  print "Using system-installed Epetra, Amesos"
#else:
import Epetra
import Amesos

def main():
  Comm = Epetra.PyComm()
  numProc = Comm.NumProc()
  iAmRoot = Comm.MyPID() == 0

  args = sys.argv[1:]
  if len(args) == 0:
    Type = "Amesos_Lapack"
  else:
    Type = args[0]

  NumGlobalRows = 10
  Map = Epetra.Map(NumGlobalRows, 0, Comm)
  LHS_exact = Epetra.MultiVector(Map, 1)
  LHS = Epetra.MultiVector(Map, 1)
  RHS = Epetra.MultiVector(Map, 1)
  Matrix = Epetra.CrsMatrix(Epetra.Copy, Map, 0)
  Indices = Epetra.IntSerialDenseVector(3)
  Values = Epetra.SerialDenseVector(3)
  Values[0] =  2.0
  Values[1] = -1.0
  Values[2] = -1.0

  NumLocalRows = Map.NumMyElements()

  # Builds the matrix (1D Laplacian)
  for ii in range(0, NumLocalRows):
    i = Map.GID(ii)
    Indices[0] = i
    if i == 0:
      NumEntries = 2
      Indices[1] = i + 1
    elif i == NumGlobalRows - 1:
      NumEntries = 2
      Indices[1] = i - 1
    else:
      NumEntries = 3
      Indices[1] = i - 1
      Indices[2] = i + 1
    Matrix.InsertGlobalValues(i, NumEntries, Values, Indices)
  ierr = Matrix.FillComplete()

  LHS_exact.Random()
  Matrix.Multiply(False, LHS_exact, RHS)
  LHS.PutScalar(1.0)

  Problem = Epetra.LinearProblem(Matrix, LHS, RHS)

  if Type == "Amesos_Lapack":
    Solver = Amesos.Lapack(Problem)
  elif Type == "Amesos_Klu":
    Solver = Amesos.Klu(Problem)
  elif Type == "Amesos_Umfpack":
    Solver = Amesos.Umfpack(Problem)
  elif Type == "Amesos_Superlu":
    Solver = Amesos.Superlu(Problem)
  elif Type == "Amesos_Superludist":
    Solver = Amesos.Superludist(Problem)
  elif Type == "Amesos_Dscpack":
    Solver = Amesos.Dscpack(Problem)
  elif Type == "Amesos_Mumps":
    Solver = Amesos.Mumps(Problem)
  else:
    print 'Selected solver (%s) not available' % Type
    print __doc__
    sys.exit(-2)
  
  AmesosList = {
    "PrintStatus": True,
    "PrintTiming": True
  }
  Solver.SetParameters(AmesosList)
  if Comm.MyPID() == 0:
    print "1) Performing symbolic factorizations..."
  Solver.SymbolicFactorization()
  if Comm.MyPID() == 0:
    print "2) Performing numeric factorizations..."
  Solver.NumericFactorization()
  if Comm.MyPID() == 0:
    print "3) Solving the linear system..."
  ierr = Solver.Solve()
  if Comm.MyPID() == 0:
    print "   Solver.Solve() return code = ", ierr
  del Solver

  # Exit with a code that indicates the total number of successes
  successes = Comm.SumAll(1)[0]
  if successes == numProc and iAmRoot: print "End Result: TEST PASSED"
  sys.exit(numProc-successes)

# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
    main()

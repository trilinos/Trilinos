#! /usr/bin/env python

# @header
# ************************************************************************
#
#              pytrilinos.amesos: python interface to amesos
#                   copyright (2005) sandia corporation
#
# under terms of contract de-ac04-94al85000, there is a non-exclusive
# license for use of this work by or on behalf of the u.s. government.
#
# this library is free software; you can redistribute it and/or modify
# it under the terms of the gnu lesser general public license as
# published by the free software foundation; either version 2.1 of the
# license, or (at your option) any later version.
#
# this library is distributed in the hope that it will be useful, but
# without any warranty; without even the implied warranty of
# merchantability or fitness for a particular purpose.  see the gnu
# lesser general public license for more details.
#
# you should have received a copy of the gnu lesser general public
# license along with this library; if not, write to the free software
# foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307
# usa
# questions? contact michael a. heroux (maherou@sandia.gov)
#
# ************************************************************************
# @header

def PrintHelp():
  print 'Usage is: ./exSolvers.py <solver-type>'
  print 'where <solver-type> can be: Amesos_Lapack, Amesos_Klu'
  print 'Amesos_Umfpack, Amesos_Superlu'

def main():
  from PyTrilinos import Amesos, Epetra

  import sys

  args = sys.argv[1:]
  if len(args) == 0:
    PrintHelp();
    sys.exit(-2);
  else:
    Type = args[0];

  n = 10;
  Comm = Epetra.SerialComm();
  Map = Epetra.Map(n, 0, Comm);
  LHS_exact = Epetra.MultiVector(Map, 1);
  LHS = Epetra.MultiVector(Map, 1);
  RHS = Epetra.MultiVector(Map, 1);
  Matrix = Epetra.CrsMatrix(Epetra.Copy, Map, 0);
  Indices = Epetra.IntSerialDenseVector(3);
  Values = Epetra.SerialDenseVector(3);
  Values[0] = 2.0; Values[1] = -1.0; Values[2] = -1.0;

  # Builds the matrix (1D Laplacian)
  for i in range(0, n):
    Indices[0] = i;
    if i == 0:
      NumEntries = 2;
      Indices[1] = i + 1;
    elif i == n - 1:
      NumEntries = 2;
      Indices[1] = i - 1;
    else:
      NumEntries = 3;
      Indices[1] = i - 1;
      Indices[2] = i + 1;
    Matrix.InsertGlobalValues(i, NumEntries, Values, Indices);
  ierr = Matrix.FillComplete();

  LHS_exact.Random();
  Matrix.Multiply(False, LHS_exact, RHS);
  LHS.PutScalar(1.0);

  Problem = Epetra.LinearProblem(Matrix, LHS, RHS);

  if Type == "Amesos_Lapack":
    Solver = Amesos.Lapack(Problem);
  elif Type == "Amesos_Klu":
    Solver = Amesos.Klu(Problem);
  elif Type == "Amesos_Umfpack":
    Solver = Amesos.Umfpack(Problem);
  elif Type == "Amesos_Superlu":
    Solver = Amesos.Superlu(Problem);
  elif Type == "Amesos_Superludist":
    Solver = Amesos.Superludist(Problem);
  elif Type == "Amesos_Dscpack":
    Solver = Amesos.Dscpack(Problem);
  elif Type == "Amesos_Mumps":
    Solver = Amesos.Mumps(Problem);
  else:
    print 'Selected solver (%s) not available' % Type
    PrintHelp();
    sys.exit(-2);
  
  AmesosList = {
    "PrintOutput": ("bool", "true"), 
    "PrintTiming": ("bool", "true"),
  };
  Solver.SetParameters(AmesosList);
  Solver.SymbolicFactorization();
  Solver.NumericFactorization();
  Solver.Solve();
  del Solver

# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
    main();

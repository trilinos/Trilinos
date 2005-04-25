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

def main(Type):
  from PyTrilinos import Amesos

  n = 10;
  Comm = Amesos.SerialComm();
  Map = Amesos.Map(n, 0, Comm);
  LHS = Amesos.MultiVector(Map, 1);
  RHS = Amesos.MultiVector(Map, 1);
  Matrix = Amesos.CrsMatrix(Amesos.Copy, Map, 0);

  # Builds the matrix (1D Laplacian)
  for i in range(0, n):
    if i != 0:
      Matrix.InsertGlobalValue(i, i - 1, -1.0);
    if i != n - 1:
      Matrix.InsertGlobalValue(i, i + 1, -1.0);
    Matrix.InsertGlobalValue(i, i, 2.0);
  ierr = Matrix.FillComplete();

  # Builds a solution that is `i' at node `i', then the
  # corresponding right-hand side, then set the solution to 0
  for i in range(0, n):
    LHS.Set(0, i, 1.0 * i);
  Matrix.Multiply(False, LHS, RHS);
  LHS.PutScalar(0.0);

  # Creates an Epetra_LinearProblem
  Problem = Amesos.LinearProblem(Matrix, LHS, RHS);

  Factory = Amesos.Factory();

  if Factory.Query(Type) == False:
    print "Selected solver (%s) not supported" % (Type)
    return;

  Solver = Factory.Create(Type, Problem);
  ierr = Solver.SetBool("PrintTiming", True);
  ierr = Solver.Solve();
  del Solver;

  error = 0.0;
  for i in range(0, n):
    error = error + abs(LHS.Get(0, i) - i);
  print "Using %s, ||x - x_ex||_1 = %e" % (Type, error);

# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
    main("Amesos_Lapack")

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

# This example creates a (distributed) tridiagonal matrix, a vector for the
# solution and one for the right-handi. Then, it solves the corresponding
# linear system using Amesos' LAPACK. Note, however, that any Amesos supported
# (and enabled at configure time) solver can be adopted; simply change the
# value of the variable `Type'
#
# Last updated on 25-Jul-5
# Author Marzio Sala, SNL 9214

# PyTrilinos imports
try:
  from PyTrilinos import Amesos, Epetra
except ImportError:
  raise ImportError, "error w/ Amesos or Epetra"

# Initializes MPI (or do-nothing in serial)
Epetra.Init()
# dimension of the problem
NumGlobalRows = 10
Comm = Epetra.PyComm()
Map = Epetra.Map(NumGlobalRows, 0, Comm)
LHS_exact = Epetra.Vector(Map)
LHS = Epetra.Vector(Map)
RHS = Epetra.Vector(Map)
Matrix = Epetra.CrsMatrix(Epetra.Copy, Map, 0)

NumLocalRows = Map.NumMyElements()

# Populates the matrix by inserting one row at-a-time. Indices and Values
# are defined as Python's lists (of the same length).
for ii in range(0, NumLocalRows):
  i = Map.GID(ii)
  if i == 0:
    Indices = [i, i + 1];
    Values  = [2.0, -1.0];
  elif i == NumGlobalRows - 1:
    Indices = [i, i - 1];
    Values  = [2.0, -1.0];
  else:
    Indices = [  i,  i - 1, i + 1];
    Values  = [2.0,   -1.0,  -1.0];
  Matrix.InsertGlobalValues(i, Values, Indices);
ierr = Matrix.FillComplete();

# Builds a solution that is `i' at node `i', then the
# corresponding right-hand side, then set the solution to 0
for i in range(0, NumLocalRows):
  LHS[i] = i;
Matrix.Multiply(False, LHS, RHS);
LHS.PutScalar(0.0);

Problem = Epetra.LinearProblem(Matrix, LHS, RHS);
Factory = Amesos.Factory();

# Creates the solver using the Amesos' factory
Type = "Amesos_Lapack"
if Factory.Query(Type) == False:
  print "Selected solver (%s) not supported" % (Type)
  Epetra.Finalize()
  raise "Solver not supported"
Solver = Factory.Create(Type, Problem);

# Setting parameters using a Python' dictionary. The list of supported
# parameters can be found on the user's guide.
AmesosList = {
  "PrintTiming":         True,
  "PrintStatus":         True,
}
Solver.SetParameters(AmesosList);

# Note: we don't check here the return parameters for brevity. 
Solver.SymbolicFactorization()
Solver.NumericFactorization()
Solver.Solve();

del Solver

error = 0.0;
for i in range(0, NumLocalRows):
  error = error + abs(LHS[i] - i);
print "Using %s, ||x - x_ex||_1 = %e" % (Type, error);

Epetra.Finalize();

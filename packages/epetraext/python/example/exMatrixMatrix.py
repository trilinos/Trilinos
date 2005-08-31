#! /usr/bin/env python

import sys

try:
  import setpath
  import Epetra, EpetraExt
except ImportError:
  from PyTrilinos import Epetra, EpetraExt
  print >>sys.stderr, "Using system-installed Epetra, EpetraExt"

def main():
  comm = Epetra.PyComm()
  n    = 10 * comm.NumProc()
  map  = Epetra.Map(n, 0, comm)

  # ================================================================= #
  # Creates two matrices, one is diagonal (A), the other contains the #
  # first sup- and super-diagonal (B), then sum them (B = A + B).     #
  # Note that B cannot be FillComplete()'d before calling Add()       #
  # unless it already contains the structure of A + B.                #
  # ================================================================= #

  Indices = Epetra.IntSerialDenseVector(2)
  Values  = Epetra.SerialDenseVector(2)

  A = Epetra.CrsMatrix(Epetra.Copy, map, 0)
  Values[0] = 2.0

  for i in range(n):
    Indices[0] = i
    A.InsertGlobalValues(i, 1, Values, Indices)
  A.FillComplete()

  B = Epetra.CrsMatrix(Epetra.Copy, map, 0)
  Values[0] = -1
  Values[1] = -1

  for i in range(n):
    if i == 0:
      NumEntries = 1
      Indices[0] = i + 1
    elif i == n - 1:
      NumEntries = 1
      Indices[0] = i - 1 
    else:
      NumEntries = 2
      Indices[0] = i - 1 
      Indices[1] = i + 1 
    B.InsertGlobalValues(i, NumEntries, Values, Indices)

  C = B

  EpetraExt.Add(A, False, 1.0, B, 1.0)
  print B


if __name__ == "__main__":
  main()

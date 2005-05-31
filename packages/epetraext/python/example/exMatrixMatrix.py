#! /usr/bin/env python
try:
  import setpath
  import Triutils, EpetraExt, Epetra
except:
  try:
    from PyTrilinos import Triutils, EpetraExt, Epetra
  except ImportError:
    raise ImportError, "error w/ Triutils or EpetraExt or Epetra"

n = 10
Comm = Epetra.SerialComm()
Map = Epetra.Map(10, 0, Comm)

# ================================================================= #
# Creates two matrices, one is diagonal (A), the other contains the #
# first sup- and super-diagonal (B), then sum them (B = A + B).     #
# Note that B cannot be FillComplete()'d before calling Add()       #
# unless it already contains the structure of A + B.                #
# ================================================================= #

Indices = Epetra.IntSerialDenseVector(2);
Values = Epetra.SerialDenseVector(2);

A = Epetra.CrsMatrix(Epetra.Copy, Map, 0);
Values[0] = 2.0;

for i in range(0, n):
  Indices[0] = i
  A.InsertGlobalValues(i, 1, Values, Indices)
A.FillComplete()

B = Epetra.CrsMatrix(Epetra.Copy, Map, 0);
Values[0] = -1;
Values[1] = -1;

for i in range(0, n):
  if i == 0:
    NumEntries = 1;
    Indices[0] = i + 1;
  elif i == n - 1:
    NumEntries = 1;
    Indices[0] = i - 1; 
  else:
    NumEntries = 2;
    Indices[0] = i - 1; 
    Indices[1] = i + 1; 
  B.InsertGlobalValues(i, NumEntries, Values, Indices)

C = B;

EpetraExt.Add(A, False, 1.0, B, 1.0)
print B

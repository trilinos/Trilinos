#! /usr/bin/env python
from PyTrilinos import Triutils, EpetraExt, Epetra

n = 10
Comm = Epetra.SerialComm()
Map = Epetra.Map(10, 0, Comm)
X = Epetra.Vector(Map)
X.Random()

# ==================================================== #
# write vector X on file "x.mm" in MatrixMarket format #
# read "map.mm", put it in Map2, then check that Map2  #
# equals Map.                                          #
# ==================================================== #

EpetraExt.BlockMapToMatrixMarketFile("map.mm", Map)

(ierr, Map2) = EpetraExt.MatrixMarketFileToBlockMap("map.mm", Comm)

if Map2.SameAs(Map):
  print "I/O for BlockMap worked!"
else:
  print "I/O for BlockMap failed!"

# ==================================================== #
# write vector X on file "x.mm" in MatrixMarket format #
# read "x.mm", put it in Y, then check that Y equals X #
# ==================================================== #

EpetraExt.MultiVectorToMatrixMarketFile("x.mm", X)

(ierr, Y) = EpetraExt.MatrixMarketFileToMultiVector("x.mm", Map)
Y.Update(1.0, X, -1.0)

if Y.Norm2()[1] == 0.0:
  print "I/O for MultiVector worked!"
else:
  print "I/O for MultiVector failed!"

# ==================================================== #
# creates a simple CrsMatrix (diagonal)                #
# write matrix A on file "A.mm" in MatrixMarket format #
# read "A.mm", put it in B, then check that B equals A #
# ==================================================== #

A = Epetra.CrsMatrix(Epetra.Copy, Map, 0);
Indices = Epetra.IntSerialDenseVector(1);
Values = Epetra.SerialDenseVector(1);

for i in range(0, n):
  Indices[0] = i
  Values[0]  = i
  A.InsertGlobalValues(i, 1, Values, Indices)
A.FillComplete()

EpetraExt.RowMatrixToMatrixMarketFile("A.mm", A)

(ierr, B) = EpetraExt.MatrixMarketFileToCrsMatrix("A.mm", Map)
EpetraExt.Add(A, False, 1.0, B, -1.0)

if B.NormInf() == 0.0:
  print "I/O for CrsMatrix worked!"
else:
  print "I/O for CrsMatrix failed!"

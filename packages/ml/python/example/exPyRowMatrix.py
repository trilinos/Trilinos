#! /usr/bin/env python
#
# Example of 1D Laplacian on one processor.
# 
# This example shows how to derive class Epetra_RowMatrix in Python.
# For the sake of simplicity, the example runs only with one processor.
# The main procedure is as follows:
# - create a Python class derived from Epetra.PyRowMatrix
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

#try:
import setpath
import Epetra
import AztecOO
#except:
#  from PyTrilinos import Epetra, AztecOO
#  print "Using installed version of Epetra, AztecOO"

class Laplace1D(Epetra.PyRowMatrix):
  def __init__(self, n, Comm):
    Epetra.PyRowMatrix.__init__(self, Comm)
    self.NumRows_ = n
    self.NumCols_ = n
    self.RowMap_  = Epetra.Map(self.NumRows_, 0, Comm)
    self.ColMap_  = Epetra.Map(self.NumCols_, 0, Comm)
    self.NumEntries_ = Epetra.IArray(1)
    self.Indices_    = Epetra.IArray(3)
    self.Values_     = Epetra.DArray(3)

  def __str__(self):
    return("Laplace1D")

  def Multiply(*args):
    self         = args[0]
    UseTranspose = args[1] 
    LHS          = args[2]
    RHS          = args[3]
    Indices    = self.Indices_
    Values     = self.Values_
    NumEntries = self.NumEntries_

    n = RHS.MyLength()
    if LHS.NumVectors() != 1:
      print "this Apply() function has been implemented for a single vector"
      return(-1)

    # I need to wrap the Values() array as DArray; if I use the bracket
    # operator the code crashes...
    LHS_V = Epetra.DArray(LHS.Values(), n)
    RHS_V = Epetra.DArray(RHS.Values(), n)

    for i in xrange(self.NumRows_):
      ierr = self.ExtractMyRowCopy(i, 5, NumEntries.Values(), Values.Values(), Indices.Values())
      total = 0.0
      for j in xrange(NumEntries[0]):
        total = total + LHS_V[Indices[j]] * Values[j]
      RHS_V[i] = total  

    return(0)

  def NumMyRowEntries(*args):
    self       = args[0]
    MyRow      = args[1]
    NumEntries = Epetra.IArray(args[2], 1)
    if MyRow == 0 | MyRow == self.NumRows_:
      NumEntries[0] = 2
    else:
      NumEntries[0] = 3
    return(0);

  # Input to this function one has an int* pointer (NumEntries),
  # a double* pointer (Values) and an int* pointer(Indices); these
  # pointers are wrapped as IArray's and DArray's
  def ExtractMyRowCopy(*args):
    self    = args[0]
    MyRow   = args[1]
    Length  = args[2]
    if (Length < 3):
      return(-1)
    NumEntries = Epetra.IArray(args[3], 1)
    Values     = Epetra.DArray(args[4], Length)
    Indices    = Epetra.IArray(args[5], Length)
    n = self.NumRows_
    if MyRow == 0:
      Indices[0] = 0; Indices[1] = 1
      Values[0] = 2.0; Values[1] = -1.0
      NumEntries[0] = 2
    elif MyRow == n - 1:
      Indices[0] = n - 1; Indices[1] = n - 2
      Values[0] = 2.0; Values[1] = -1.0
      NumEntries[0] = 2
    else:
      Indices[0] = MyRow; Indices[1] = MyRow - 1; Indices[2] = MyRow + 1
      Values[0] = 2.0; Values[1] = -1.0; Values[2] = -1.0
      NumEntries[0] = 3
    return(0)

  def ApplyInverse(*args):
    return(-2)
  def HasNormInf(*args):
    return(True);

  def NormInf(*args):
    return(0)
  def MaxNumEntries(*args):
    return(3)
  def Apply(*args):
    self = args[0]
    LHS = args[1]
    RHS = args[2]
    return(self.Multiply(False, LHS, RHS))
  def NumGlobalNonzeros(*args):
    self = args[0]
    return(3 * self.NumRows_ - 2);
  def NumMyNonzeros(*args):
    self = args[0]
    return(3 * self.NumRows_ - 2);
  def NumGlobalDiagonals(*args):
    self = args[0]
    return(self.NumRows_);
  def NumMyDiagonals(*args):
    self = args[0]
    return(self.NumRows_);
  def ExtractDiagonalCopy(*args):
    self = args[0]
    Diagonal = args[1]
    Diagonal.PutScalar(2.0)
    return(0)
  def HasNormInf(*args):
    return(True)
  def NormInf(*args):
    return(4.0)
  def NormOne(*args):
    return(4.0)
  def UseTranspose(*args):
    return(False)
  def SetTranspose(*args):
    return(0)

  def RowMatrixRowMap(*args):
    self = args[0]
    return(self.RowMap_);

  def RowMatrixColMap(*args):
    self = args[0]
    return(self.ColMap_);

  def UpperTriangular(*args):
    return(False)

  def LowerTriangular(*args):
    return(False)

  def NumMyRows(*args):
    self = args[0]
    return(self.NumRows_);

  def NumMyCols(*args):
    self = args[0]
    return(self.NumCols_);

  def NumGlobalRows(*args):
    self = args[0]
    return(self.NumRows_);

  def NumGlobalCols(*args):
    self = args[0]
    return(self.NumCols_);

  def Filled(*args):
    return(True)

  def RightScale(*args):
    return(-1)

  def LeftScale(*args):
    return(-1)

  def InvRowSum(*args):
    return(-1)

  def InvColSum(*args):
    return(-1)

  def Solve(*args):
    return(-1)

  def OperatorDomainMap(*args):
    self = args[0]
    return(self.RowMap_)

  def OperatorRangeMap(*args):
    self = args[0]
    return(self.RowMap_)

  def Map(*args):
    self = args[0]
    return(self.RowMap_)

  def RowMatrixImporter(*args):
    self = args[0]
    Importer = Epetra.Import(self.RowMap_, self.RowMap_)
    return(Importer)

def main():

  n = 10
  Comm = Epetra.PyComm()
  if Comm.NumProc() != 1:
    print "This example is only serial, sorry"
    return
  Matrix = Laplace1D(n, Comm)

  LHS = Epetra.Vector(Matrix.Map())
  RHS = Epetra.Vector(Matrix.Map())
  LHS.PutScalar(0.0)
  RHS.PutScalar(1.0)
  
  Solver = AztecOO.AztecOO(Matrix, LHS, RHS)
  
  Solver.SetAztecOption(AztecOO.AZ_solver,AztecOO.AZ_gmres)
  Solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_Jacobi)
  Solver.SetAztecOption(AztecOO.AZ_output, 16)
  Solver.Iterate(1550, 1e-5)
  
  del Matrix

# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX
# command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()

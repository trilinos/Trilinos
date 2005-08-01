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
#       pointers or references. Some tools are available to do that;
#       for example a new int array is allocated using NewInt(),
#       an element is set using SetInt() and get using GetInt(),
#       and the array is destroy'd using DeleteInt() -- equivalently
#       for double arrays.
#
# \author Marzio Sala, 9214
#
# \date Last updated on 31-Jul-05

try:
  import setpath
  import Epetra, AztecOO
except:
  try:
    from PyTrilinos import Epetra, AztecOO
  except ImportError:
    raise ImportError, "error w/ Epetra or AztecOO"

class Laplace1D(Epetra.PyRowMatrix):
  def __init__(self, n, Comm):
    Epetra.PyRowMatrix.__init__(self, Comm)
    self.NumRows_ = n
    self.NumCols_ = n
    self.RowMap_  = Epetra.Map(self.NumRows_, 0, Comm)
    self.ColMap_  = Epetra.Map(self.NumCols_, 0, Comm)
    self.NumEntries_ = Epetra.NewInt(1)
    self.Indices_    = Epetra.NewInt(3)
    self.Values_     = Epetra.NewDouble(3)

  def __del__(self):
    Epetra.DeleteInt(self.NumEntries_);
    Epetra.DeleteInt(self.Indices_);
    Epetra.DeleteDouble(self.Values_);
   
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

    for i in xrange(self.NumRows_):
      ierr = self.ExtractMyRowCopy(i, 5, NumEntries, Values, Indices)
      total = 0.0
      for j in xrange(Epetra.GetInt(NumEntries,0)):
        total = total + LHS[0,Epetra.GetInt(Indices,j)] * Epetra.GetDouble(Values,j)
      RHS[0,i] = total  

    return(0)

  def NumMyRowEntries(*args):
    self       = args[0]
    MyRow      = args[1]
    NumEntries = args[2]
    if MyRow == 0 | MyRow == self.NumRows_:
      Epetra.SetInt(NumEntries, 0, 2)
    else:
      Epetra.SetInt(NumEntries, 0, 3)
    return(0);

  def ExtractMyRowCopy(*args):
    self    = args[0]
    MyRow   = args[1]
    Length  = args[2]
    if (Length < 3):
      return(-1)
    NumEntries = args[3]
    Values     = args[4]
    Indices    = args[5]
    n = self.NumRows_
    if MyRow == 0:
      Epetra.SetDouble(Values, 0, 2.0)
      Epetra.SetDouble(Values, 1, -1.0)
      Epetra.SetInt(Indices, 0, 0)
      Epetra.SetInt(Indices, 1, 1)
      Epetra.SetInt(NumEntries, 0, 2)
    elif MyRow == n - 1:
      Epetra.SetDouble(Values, 0, 2.0)
      Epetra.SetDouble(Values, 1, 1.0)
      Epetra.SetInt(Indices, 0, n - 1)
      Epetra.SetInt(Indices, 1, n - 2)
      Epetra.SetInt(NumEntries, 0, 2)
    else:
      Epetra.SetDouble(Values, 0, 2.0)
      Epetra.SetDouble(Values, 1, -1.0)
      Epetra.SetDouble(Values, 2, -1.0)
      Epetra.SetInt(Indices, 0, MyRow)
      Epetra.SetInt(Indices, 1, MyRow - 1)
      Epetra.SetInt(Indices, 2, MyRow + 1);
      Epetra.SetInt(NumEntries, 0, 3)
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

  def Importer(*args):
    return(Null)

def main():

  n = 10
  Comm = Epetra.PyComm()
  if Comm.NumProc() != 1:
    print "This example is only serial, sorry"
    return
  Matrix = Laplace1D(n, Comm)

  # Use Multivectors, not vectors!
  LHS = Epetra.MultiVector(Matrix.Map(), 1)
  RHS = Epetra.MultiVector(Matrix.Map(), 1)
  LHS. PutScalar(0.0)
  RHS. PutScalar(1.0)
  
  Problem = Epetra.LinearProblem(Matrix, LHS, RHS)
  Solver = AztecOO.AztecOO(Problem)
  
  Solver.SetAztecOption(AztecOO.AZ_solver,AztecOO.AZ_gmres)
  Solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_Jacobi)
  Solver.SetAztecOption(AztecOO.AZ_output, 16)
  #Solver.SetPrecOperator(Prec)
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

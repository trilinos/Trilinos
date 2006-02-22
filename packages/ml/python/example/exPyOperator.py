#! /usr/bin/env python
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

try:
  import setpath
  import Epetra
  import AztecOO
except:
  from PyTrilinos import Epetra, AztecOO
  print "Using installed versions of Epetra, AztecOO"

############################################################

class MyOperator(Epetra.PyOperator):

  def __init__(self, map):
    Epetra.PyOperator.__init__(self, map.Comm())
    self.__map   = map
    self.__label = "1D Laplace"

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
    n   = rhs.MyLength()
    rhs[ 0 ] = 2.0 * lhs[ 0 ] - lhs[ 1 ]
    rhs[n-1] = 2.0 * lhs[n-1] - lhs[n-2]
    for i in xrange(1, n-1):
      rhs[i] = 2.0 * lhs[i] - lhs[i-1] - lhs[i+1]
    return 0

  def ApplyInverse(self):
    return -1

  def OperatorDomainMap(self):
    return self.__map

  def OperatorRangeMap(self):
    return self.__map

  def Map(self):
    return self.__map

  def HasNormInf(self):
    return True

  def NormInf(self):
    return 0

############################################################

def main():

  n    = 100
  comm = Epetra.Pycomm()
  if comm.NumProc() != 1:
    print "This example is only serial, sorry"
    return
  map = Epetra.Map(n, 0, comm)
  op  = MyOperator(map)

  print op.Label()
  lhs = Epetra.Vector(map)
  rhs = Epetra.Vector(map)
  
  rhs.PutScalar(1.0)
  lhs.PutScalar(0.0)
  
  Problem = Epetra.LinearProblem(op, lhs, rhs)
  Solver  = AztecOO.AztecOO(Problem)
  
  Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg)
  Solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_none)
  Solver.SetAztecOption(AztecOO.AZ_output, 16)
  Solver.Iterate(1550, 1e-5)
  
# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX
# command line.  
# This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()


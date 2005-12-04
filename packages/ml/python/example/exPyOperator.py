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
# \author Marzio Sala, 9214
#
# \date Last updated on 31-Jul-05

try:
  import setpath
  import Epetra
  import AztecOO
except:
  from PyTrilinos import Epetra, AztecOO
  print "Using installed versions of Epetra, AztecOO"

class MyOperator(Epetra.PyOperator):
  def __init__(self, Map):
    Epetra.PyOperator.__init__(self, Map.Comm());
    self.Map_ = Map

  def __str__(self):
    return("MyPyOperator")

  def Apply(*args):
    LHS = args[1]
    RHS = args[2]
    n = RHS.MyLength()
    if LHS.NumVectors() != 1:
      print "this Apply() function has been implemented for a single vector"
      return(-1)

    LHS_V = Epetra.DArray(LHS.Values(), n)
    RHS_V = Epetra.DArray(RHS.Values(), n)

    # I don't undestand why, but if I use the bracket operator the
    # code crashes...
    RHS_V[0] = 2.0 * LHS_V[0] - LHS_V[1]
    RHS_V[n - 1] = 2.0 * LHS_V[n - 1] - LHS_V[n - 2]
    for i in xrange(1, n - 1):
      RHS_V[i] = 2.0 * LHS_V[i] - LHS_V[i - 1] - LHS_V[i + 1]
    return(0)

  def ApplyInverse(*args):
    return(-1)

  def OperatorDomainMap(*args):
    self = args[0]
    return(self.Map_)

  def OperatorRangeMap(*args):
    self = args[0]
    return(self.Map_)

  def Map(*args):
    self = args[0]
    return(self.Map_)

  def HasNormInf(*args):
    return(True);

  def NormInf(*args):
    return(0)

def main():

  n = 100
  Comm = Epetra.PyComm()
  if Comm.NumProc() != 1:
    print "This example is only serial, sorry"
    return
  Map = Epetra.Map(n, 0, Comm)
  Op = MyOperator(Map)

  print Op.Label()
  LHS = Epetra.Vector(Map)
  RHS = Epetra.Vector(Map)
  
  RHS.PutScalar(1.0);
  LHS.PutScalar(0.0);
  
  Problem = Epetra.LinearProblem(Op, LHS, RHS)
  Solver = AztecOO.AztecOO(Problem)
  
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


#! /usr/bin/env python

# --------------------------------------------------------------------------- #
# Simple class that defines a multilevel preconditioner based on aggregation.
# The class requires in input the matrix (serial or distributed) defined as an
# ML.Operator, and the maximum number of levels. 
#
# NOTE: This code does not check for the actual number of levels; besides,
#       the smoother is always symmetric Gauss-Seidel. The nullspace is
#       supposed to be constant, with only one unknown. Simple changes can
#       be made to increase the flexibility of the code. If you want to 
#       define your own smoother, check example exMLAPI_Smoother.py
#
# \author Marzio Sala, SNL 9214
#
# \date Last updated on 05-Aug-05
# --------------------------------------------------------------------------- #

from   optparse import *
import sys

parser = OptionParser()
parser.add_option("-t", "--testharness", action="store_true",
                  dest="testharness", default=False,
                  help="test local build modules; prevent loading system-installed modules")
parser.add_option("-v", "--verbosity", type="int", dest="verbosity", default=2,
                  help="set the verbosity level [default 2]")
options,args = parser.parse_args()
if options.testharness:
  import setpath
  import Epetra
  import ML
else:
  try:
    import setpath
    import Epetra
    import ML
  except ImportError:
    from PyTrilinos import Epetra
    from PyTrilinos import ML
    print >>sys.stderr, "Using system-installed Epetra, ML"

################################################################################

class MultiLevel:
  def __init__(self, A, MaxLevels):
    
    # Declares the lists that will contain the matrices (A), the prolongator
    # operators (P), the restriction operator (R), and the smoother or
    # coarse solver (S). The finest level matrix will be stored in A[0].
    A_array = [A]  
    P_array = []
    R_array = []
    S_array = []

    self.MaxLevels_ = MaxLevels;
    self.NullSpace_ = ML.MultiVector(A.GetDomainSpace())
    self.NullSpace_.Update(1.0)

    for level in xrange(MaxLevels):
      # Needs a vector containing the actual null space, plus another
      # vector for the next-level null space.
      ThisNS = self.NullSpace_
      NextNS = ML.MultiVector()
      A = A_array[level]
      # Constructs the non-smoothed prolongator...
      List = {
        "aggregation: type": "MIS",
        "aggregation: threshold": 0.00,
        "PDE equations": 1
      }
      Ptent = ML.GetPNonSmoothed(A, ThisNS, NextNS, List)
      # ...and smooth it. Here we skip the diagonal scaling
      I = ML.GetIdentity(A.GetDomainSpace(), A.GetRangeSpace())
      lambda_max = 1.0 / ML.MaxEigAnorm(A)
      P = (I - A * lambda_max) * Ptent
      # Stores prolongator, restriction, and RAP product
      P = (I - A * lambda_max) * Ptent
      R = ML.GetTranspose(P)
      A_coarse = ML.GetRAP(R, A, P)
      # Defines the coarse solver or smoother (symmetric Gauss-Seidel)
      if level == MaxLevels - 1:
        S = ML.InverseOperator(A, "Amesos")
      else:
        S = ML.InverseOperator(A, "symmetric Gauss-Seidel")
      
      A_array.append(A_coarse)
      P_array.append(P)
      R_array.append(R)
      S_array.append(S)
      self.NullSpace_ = NextNS
    
    self.A_array_ = A_array
    self.P_array_ = P_array
    self.R_array_ = R_array
    self.S_array_ = S_array

  def Apply(self, RHS):
    LHS = self.MultiLevelCycle(RHS, 0)
    return(LHS)
         
  def MultiLevelCycle(self, b_f, level):
    A = self.A_array_[level];
    P = self.P_array_[level];
    R = self.R_array_[level];
    S = self.S_array_[level];
    MaxLevels = self.MaxLevels_

    if level == MaxLevels - 1:
      return(S * b_f)

    # apply pre-smoother
    x_f = S * b_f
    # new residual
    r_f = b_f - A * x_f
    # restrict to coarse
    r_c = R * r_f
    # solve coarse problem
    z_c = self.MultiLevelCycle(r_c, level + 1)
    # prolongate back and add to solution
    x_f = x_f + P * z_c
    # apply post-smoother
    S.Apply(b_f, x_f)
  
    return(x_f)

################################################################################

def main():
  # Defines a communicator (serial or parallel, depending on how Trilinos
  # was configured), and creates a matrix corresponding to a 1D Laplacian.
  # AT THIS MOMENT THE EXAMPLE IS ONLY SERIAL
  comm  = Epetra.PyComm()
  if comm.NumProc() > 1: return

  n     = 1000
  Space = ML.Space(n)

  Matrix = ML.PyMatrix(Space, Space)

  for i in Space.GetMyGlobalElements():
    if i > 0:
      Matrix[i, i - 1] = -1.
    if i < n - 1:
      Matrix[i, i + 1] = -1.
    Matrix[i, i] = 2.0
  
  Matrix.FillComplete()
  
  MaxLevels = 3
  Prec = MultiLevel(Matrix, MaxLevels)
  
  # Defines a linear system with a given exact solution
  x = ML.MultiVector(Space)
  x_exact = ML.MultiVector(Space)
  x_exact.Random()
  x_exact_norm = x_exact.Norm2()
  y = Matrix * x_exact
  
  # ---------------------------------------------------- #
  # A very simple Richardson method to solve the problem #
  # using 10 iterations:                                 #
  # ---------------------------------------------------- #
  for i in xrange(10):
    # compute non-preconditioned residual
    res = y - Matrix * x 
    # compute the preconditioned residual
    prec_res = Prec.Apply(res)
    # update the solution vector
    x = x + prec_res 
    # compute the energy of the error
    diff = x - x_exact
    norm = diff * diff  # or use diff * (Matrix * diff)
    if comm.MyPID() == 0:
      print "iter ", i, " ||x - x_exact||_2 = ", norm
      print "End Result: TEST PASSED"

# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()

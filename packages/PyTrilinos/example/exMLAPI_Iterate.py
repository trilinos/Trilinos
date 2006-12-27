#! /usr/bin/env python

# @HEADER
# ************************************************************************
# 
#                PyTrilinos: Python Interface to Trilinos
#                   Copyright (2005) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Bill Spotz (wfspotz@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# --------------------------------------------------------------------------- #
# Simple class that defines a multilevel preconditioner based on aggregation.
# The class requires in input the matrix (serial or distributed) defined as an
# ML.Operator, and the maximum number of levels. 
#
# NOTE: This code does not check for the actual number of levels; besides,
#       the smoother is always symmetric Gauss-Seidel. Simple changes can
#       be made to increase the flexibility of the code. If you want to 
#       define your own smoother, check example exMLAPI_Smoother.py
#
# \author Marzio Sala, SNL 9214
#
# \date Last updated on 03-Aug-05
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
  except:
    from PyTrilinos import Epetra
    from PyTrilinos import ML
    print >>sys.stderr, "Using system-installed versions of Epetra, ML"

################################################################################

class MultiLevel(ML.BaseOperator):
  def Reshape(self, Matrix, MaxLevels):
    
    self.SetLabel("MultiLevel Preconditioner")
    # Declares the lists that will contain the matrices (A), the prolongator
    # operators (P), the restriction operator (R), and the smoother or
    # coarse solver (S). The finest level matrix will be stored in A[0].
    A_array = [Matrix]  
    P_array = []
    R_array = []
    S_array = []

    self.MaxLevels_ = MaxLevels;
    self.Matrix_ = Matrix
    self.DomainSpace_ = Matrix.GetDomainSpace()
    self.RangeSpace_ = Matrix.GetRangeSpace()

    NullSpace = ML.MultiVector(Matrix.GetDomainSpace())
    NullSpace.Update(1.0)

    List = {
      "aggregation: type": "MIS"
    }

    for level in xrange(MaxLevels):
      NextNullSpace = ML.MultiVector()
      A = A_array[level]
      # Constructs the non-smoothed prolongator...
      Ptent = ML.GetPNonSmoothed(A, NullSpace, NextNullSpace, List)
      # ...and smooth it with (I - D A / lambda_max{D A}),
      # where D contains the inverse of the diagonal of A,
      # and lambda_max the maximum eigenvalue. This is done is several
      # steps: first, we extract the diagonal of A[level] as a vector,
      # then we compute the inverse of it, and we form a matrix with this
      # diagonal entries. Then, we need to allocate an identity matrix,
      # and to compute lambda_max.
      D_vector = ML.GetDiagonal(A)
      D_vector.Reciprocal()
      D = ML.GetDiagonal(D_vector)
      I = ML.GetIdentity(A.GetDomainSpace(), A.GetRangeSpace())
      lambda_max = 1.0 / ML.MaxEigAnorm(D * A)
      P = (I - D * A * lambda_max) * Ptent
      R = ML.GetTranspose(P)
      A_coarse = ML.GetRAP(R, A, P)
      # Defines the coarse solver or smoother
      if level == MaxLevels - 1:
        S = ML.InverseOperator(A, "Amesos")
      else:
        List = {
          "smoother: sweeps": 2,
          "smoother: damping factor": 0.67
        }
        S = ML.InverseOperator()
        S.Reshape(A, "Jacobi", List)
        # Other choices for smoothers are:
        # - "Gauss-Seidel"
        # - "symmetric Gauss-Seidel"
        # - "ILU", "IC", "ILUT", "ICT"

      # Stores prolongator, restriction, RAP product and smoother (or
      # coarse solver)
      P_array.append(P)
      R_array.append(R)
      A_array.append(A_coarse)
      S_array.append(S)

      NullSpace = NextNullSpace
    
    self.A_array_ = A_array
    self.P_array_ = P_array
    self.R_array_ = R_array
    self.S_array_ = S_array

  def Apply(*args):
    self = args[0]
    RHS = args[1]
    LHS = args[2]
    LHS.Update(self.MultiLevelCycle(RHS, 0))
    return(0)
         
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

  def GetOperatorDomainSpace(self):
    return(self.DomainSpace_)

  def GetOperatorRangeSpace(self):
    return(self.RangeSpace_)

  def __str__(self):
    return "MultiLevel"

  def __mul__(*args):
    self = args[0]
    rhs = args[1]
    res = ML.MultiVector(rhs.GetVectorSpace())
    self.Apply(rhs, res)
    return(res)

################################################################################

# -------------------------------------------------------------------------- #
# Main driver.
# -------------------------------------------------------------------------- #
#
def main():

  # Creates a matrix corresponding to a 1D Laplacian.
  comm = Epetra.PyComm()
  if comm.NumProc() > 1: return

  n = 1000
  FineSpace = ML.Space(n)
  
  Matrix = ML.PyMatrix(FineSpace, FineSpace)
  for i in xrange(n):
    if i > 0:
      Matrix[i, i - 1] = -1.
    if i < n - 1:
      Matrix[i, i + 1] = -1.
    Matrix[i, i] = 2.0
  
  Matrix.FillComplete()
  
  # Allocates the preconditioner
  MaxLevels = 2
  Prec = MultiLevel()
  Prec.Reshape(Matrix, MaxLevels)
  
  # Defines a linear system, and solve it using AztecOO.
  LHS = ML.MultiVector(FineSpace)
  RHS = ML.MultiVector(FineSpace)
  RHS.Random()
 
  # Iterate using AztecOO's GMRES
  List = {
    "krylov: type": "gmres",
    "krylov: tolerance": 10e-5
  }
  ML.Iterate(Matrix, LHS, RHS, Prec, List)

  if comm.MyPID() == 0: print "End Result: TEST PASSED"

################################################################################

if __name__ == "__main__":
  main()

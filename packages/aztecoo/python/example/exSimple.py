#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#             PyTrilinos.AztecOO: Python Interface to AztecOO
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
# @HEADER

from optparse import *

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
   import Galeri
   import AztecOO
else:
   try:
      import setpath
      import Epetra
      import Galeri
      import AztecOO
   except ImportError:
      from PyTrilinos import Epetra
      from PyTrilinos import Galeri
      from PyTrilinos import AztecOO
      print "Using system-installed Epetra, Galeri, AztecOO"

comm = Epetra.PyComm()

nx = 30; ny = 30
GaleriList = {"n": nx * ny,  # for Linear map
              "nx": nx,      # for Laplace2D, which requires nx
              "ny": ny       # and ny
              }
Map = Galeri.CreateMap("Linear", comm, GaleriList)
Matrix = Galeri.CreateCrsMatrix("Laplace2D", Map, GaleriList)
Exact = Epetra.Vector(Map) 
LHS = Epetra.Vector(Map) 
RHS = Epetra.Vector(Map) 
Exact.Random()       # fix exact solution
LHS.PutScalar(0.0)   # fix starting solution
Matrix.Multiply(False, Exact, RHS) # fix rhs corresponding to Exact

# Solve the linear problem
if 0:
   # this does not work on most installations
   Problem = Epetra.LinearProblem(Matrix, LHS, RHS)
   Solver = AztecOO.AztecOO(Problem)
else:
   Solver = AztecOO.AztecOO(Matrix, LHS, RHS)

try:
   plist = {"Solver"          : "CG",
            "Precond"         : "Dom_Decomp",
            "Subdomain_Solve" : "ICC",
            "Output"          : 16
            }
   Solver.SetParameters(plist,True)
except AttributeError:
   # If AztecOO and its python wrappers have been compiled without Teuchos
   # support, then then the SetParameters method will not exist, thus raising an
   # AttributeError exception
   print "Teuchos ParameterLists not supported"
   Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg)
   Solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_dom_decomp)
   Solver.SetAztecOption(AztecOO.AZ_subdomain_solve, AztecOO.AZ_icc)
   Solver.SetAztecOption(AztecOO.AZ_output, 16)
   Solver.Iterate(1550, 1e-5)

if comm.MyPID() == 0: print "End Result: TEST PASSED"

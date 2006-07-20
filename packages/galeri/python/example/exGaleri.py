#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#              PyTrilinos.Galeri: Python Interface to Galeri
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

# "from PyTrilinos import ..." syntax.  Here, the setpath module adds the build
# directory, including "PyTrilinos", to the front of the search path.  We thus
# use "import ..." for Trilinos modules.  This prevents us from accidentally
# picking up a system-installed version and ensures that we are testing the
# build module.
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
    import Galeri
else:
    try:
        import setpath
        import Epetra
        import Galeri
    except ImportError:
        from PyTrilinos import Epetra, Galeri
        print >>sys.stderr, "Using system-installed Epetra, Galeri"

# Create a communicator
comm    = Epetra.PyComm()
numProc = comm.NumProc()
iAmRoot = comm.MyPID() == 0

# Create a Cartesian map, containing nx x ny x NumProcs nodes
nx = 2
ny = 2 * numProc
List = {"nx": nx,               # number of nodes in the X-direction
        "ny": ny,               # number of nodes in the Y-directioN
        "mx": 1,                # number of processors in the X-direction
        "my": numProc           # number of processors in the Y-direction
        }

# Creating a map corresponding to a 2D Cartesian distribuion
if iAmRoot: print "Creating a Map..."

Map = Galeri.CreateMap("Cartesian2D", comm, List)

# creates a point matrix based on the previously created map
if iAmRoot: print "Creating an Epetra_CrsMatrix..."

CrsMatrix = Galeri.CreateCrsMatrix("Laplace2D", Map, List)

# extend the matrix into VBR format, by replicating each equation (twice
# in this example)
if iAmRoot: print "Extending the Epetra_CrsMatrix into an Epetra_VbrMatrix..."

VbrMatrix = Galeri.CreateVbrMatrix(CrsMatrix, 2);

# Now work a bit with vectors
LHS = Epetra.Vector(Map)
LHS.PutScalar(12.0)
RHS = Epetra.Vector(Map)
#LinearProblem = Epetra.LinearProblem(CrsMatrix, LHS, RHS)
# Now the linear problem is not solved, use for example AztecOO with IFPACK or
# ML, or Amesos
norm = Galeri.ComputeNorm(CrsMatrix, LHS, RHS)
if iAmRoot: print "||A x - b||_2 = ", norm

# Indicate success
successes = comm.SumAll(1)
if successes == numProc and iAmRoot: print "End Result: TEST PASSED"
sys.exit(numProc-successes)

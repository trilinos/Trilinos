#! /usr/bin/env python

# "from PyTrilinos import ..." syntax.  Here, the setpath module adds the build
# directory, including "PyTrilinos", to the front of the search path.  We thus
# use "import ..." for Trilinos modules.  This prevents us from accidentally
# picking up a system-installed version and ensures that we are testing the
# build module.
try:
  import setpath
  import Epetra, Galeri
except ImportError:
  from PyTrilinos import Epetra, Galeri
  print "Using system-installed Epetra, Galeri"

Comm = Epetra.PyComm()
# Create a Cartesian map, containing nx x ny x NumProcs nodes
nx = 2
ny = 2 * Comm.NumProc()
List = {
  "nx": nx,               # number of nodes in the X-direction
  "ny": ny,               # number of nodes in the Y-directioN
  "mx": 1,                # number of processors in the X-direction
  "my": Comm.NumProc()    # number of processors in the Y-direction
}

# Creating a map corresponding to a 2D Cartesian distribuion
if Comm.MyPID() == 0:
  print "Creating a Map..."

Map = Galeri.CreateMap("Cartesian2D", Comm, List)

# creates a point matrix based on the previously created map
if Comm.MyPID() == 0:
  print "Creating an Epetra_CrsMatrix..."

CrsMatrix = Galeri.CreateCrsMatrix("Laplace2D", Map, List)

# extend the matrix into VBR format, by replicating each equation (twice
# in this example)
if Comm.MyPID() == 0:
  print "Extending the Epetra_CrsMatrix into an Epetra_VbrMatrix..."

VbrMatrix = Galeri.CreateVbrMatrix(CrsMatrix, 2);

# Now work a bit with vectors
LHS = Epetra.Vector(Map)
LHS.PutScalar(12.0)
RHS = Epetra.Vector(Map)
LinearProblem = Epetra.LinearProblem(CrsMatrix, LHS, RHS)
# Now the linear problem is not solved, use for example AztecOO with IFPACK or
# ML, or Amesos
norm = Galeri.ComputeNorm(CrsMatrix, LHS, RHS)
if Comm.MyPID() == 0:
  print "||A x - b||_2 = ", norm



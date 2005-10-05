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
  "nx": nx,
  "ny": ny
}

Map = Galeri.CreateMap("Cartesian2D", Comm, List)
print Map
Matrix = Galeri.CreateCrsMatrix("Laplace2D", Map, List)
print Matrix

LHS = Epetra.Vector(Map)
LHS.PutScalar(12.0)
RHS = Epetra.Vector(Map)
LinearProblem = Epetra.LinearProblem(Matrix, LHS, RHS)
# Now the linear problem is not solved, use for example AztecOO with IFPACK or
# ML, or Amesos
print Galeri.ComputeNorm(Matrix, LHS, RHS)



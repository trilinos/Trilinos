from PyTrilinos import Epetra, Galeri

Comm = Epetra.PyComm()
# Create a Cartesian map, containing nx x ny x NumProcs nodes
nx = 2
ny = 2 * Comm.NumProc()
List = {
  "n": nx * ny,
  "nx": nx,
  "ny": ny
}

Map = Galeri.CreateMap("linear", Comm, List)
print Map
Matrix = Galeri.CreateCrsMatrix("Laplace2D", Map, List)
print Matrix

LHS = Epetra.Vector(Map)
LHS.PutScalar(12.0)
RHS = Epetra.Vector(Map)
LinearProblem = Epetra.LinearProblem(Matrix, LHS, RHS)

print Galeri.ComputeNorm("2-norm", Matrix, LHS, RHS)



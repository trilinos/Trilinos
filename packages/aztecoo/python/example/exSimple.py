#! /usr/bin/env python
try:
   import setpath
except ImportError:
   from PyTrilinos import Epetra, Triutils, AztecOO
   print "Using system-installed Epetra, Triutils, AztecOO"
else:
   import Epetra
   import Triutils
   import AztecOO

nx = 100
ny = 100
Comm = Epetra.PyComm()
Gallery = Triutils.CrsMatrixGallery("laplace_2d", Comm)
Gallery.Set("nx", nx)
Gallery.Set("ny", ny)
Matrix = Gallery.GetMatrix()
LHS = Gallery.GetStartingSolution()
RHS = Gallery.GetRHS()

# Solve the linear problem
if 0:
  Problem = Epetra.LinearProblem(Matrix, LHS, RHS)
  Solver = AztecOO.AztecOO(Problem)
else:
  Solver = AztecOO.AztecOO(Matrix, LHS, RHS)

Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg)
Solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_dom_decomp)
Solver.SetAztecOption(AztecOO.AZ_subdomain_solve, AztecOO.AZ_icc)
Solver.SetAztecOption(AztecOO.AZ_output, 16)
Solver.Iterate(1550, 1e-5)

if Comm.MyPID() == 0: print "End Result: TEST PASSED"

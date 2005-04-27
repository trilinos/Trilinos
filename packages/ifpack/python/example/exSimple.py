#! /usr/bin/env python
from PyTrilinos import Triutils, IFPACK, AztecOO, Epetra
nx = 20;
ny = 20;
nz = 20;
Comm = Epetra.SerialComm();
Gallery = Triutils.CrsMatrixGallery("laplace_3d", Comm)
Gallery.Set("nx", nx);
Gallery.Set("ny", ny);
Gallery.Set("nz", nz);
Matrix = Gallery.GetMatrix();
LHS = Gallery.GetStartingSolution();
RHS = Gallery.GetRHS();

Factory = IFPACK.Factory();
Prec = Factory.Create("IC", Matrix);
Prec.SetInt("fact: level-of-fill", 5);
Prec.Initialize();
Prec.Compute();

Problem = Epetra.LinearProblem(Matrix, LHS, RHS);
Solver = AztecOO.AztecOO(Problem);

Solver.SetPrecOperator(Prec)
Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg);
Solver.SetAztecOption(AztecOO.AZ_output, 16);
Solver.Iterate(1550, 1e-5)

print Prec

#! /usr/bin/env python
from PyTrilinos import ML, Triutils, IFPACK, AztecOO, Epetra
nx = 100;
ny = 100;
Comm = Epetra.SerialComm();
Gallery = Triutils.CrsMatrixGallery("laplace_2d", Comm)
Gallery.Set("nx", nx);
Gallery.Set("ny", ny);
Matrix = Gallery.GetMatrix();
LHS = Gallery.GetStartingSolution();
RHS = Gallery.GetRHS();

Factory = IFPACK.Factory();
Prec = ML.MultiLevelPreconditioner(Matrix, False);
Prec.SetInt("max levels", 3);
Prec.SetString("smoother: type (level 0)", "Aztec");
Prec.ComputePreconditioner();
Solver = AztecOO.AztecOO(Matrix, LHS, RHS)

Solver.SetPrecOperator(Prec)
Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg);
Solver.SetAztecOption(AztecOO.AZ_output, 16);
Solver.Iterate(1550, 1e-5)

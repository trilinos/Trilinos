To Compile and HyperLU
------------------

1. Compile with Trilinos with the following packages enabled. Zoltan,
Isorropia, Epetra, EpetraExt, AztecOO, Amesos, Ifpack, Teuchos and Belos.

Another option is to compile Trilinos with TrilinosCouplings with Isorropia
and Belos. This will compile more packages than required by HyperLU, but it is
an easy way to get started.

2. Copy the file Trilinos/BUILD_DIR/Makefile.export.Trilinos_install to
ILU/epetra.

3. Run "make" or "make all" to compile the drivers for HyperLU. There is a
AztecOO based driver (hyperlu_driver.exe) and a Belos based driver
(hyperlu_belos_driver.exe).

4. Edit HyperLU.xml to set the precondtioner (HyperLU)
and the matrix market file name.
    Please leave the parameters "Schur Approximation Method" and "Symmetry" as
    it is. I will end up doing a symmetric permutation of an unsymmetric
    matrix for now. The approximation method is currently A22 + a percentage
    of the block diagonals. Adjust the diagonal factor (0.0 to 1.0) for 
    better convergence (at the cost of slower setup and more memory usage)

    The inner solver is currently AztecOO with a domain decomposition 
    preconditioner. It is not expose in the interface yet. I can do that if
    you need to change that.

5. mpirun -n np <driver_name> will solve the linear system from the
matrix market file with a right hand side that is all ones. This is probably
the best way to run for now if you have the matrices in a matrix market file. 
Please use Ifpack_HyperLU as your preconditioner otherwise. 

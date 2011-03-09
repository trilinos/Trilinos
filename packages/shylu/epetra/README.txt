To Compile and HyperLU
------------------

1. Compile with Trilinos with TrilinosCouplings and Isorropia enabled. This will
compile more packages than required by HyperLU, but it is an easy way to get
started.

2. Copy the file Trilinos/BUILD_DIR/Makefile.export.Trilinos_install to
ILU/epetra.

3. Run "make" to compiler the drivers for HyperLU. There is a AztecOO based
driver (hyperlu_driver.exe) and a Belos based driver
(hyperlu_belos_driver.exe).

4. Edit HyperLU.xml to set the precondtioner (HyperLU/ILU/ILUT/ML currently)
and the matrix market file name.

5. mpirun -n np <driver_name> will solve the linear system from the
matrix market file with a right hand side that is all ones.

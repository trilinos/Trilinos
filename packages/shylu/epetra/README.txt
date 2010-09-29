To Compile and HyperLU
------------------

1. Compile with Trilinos with TrilinosCouplings and Isorropia enabled. This will
compile more packages than required by HyperLU, but it is an easy way to get
started.

2. Copy the file Trilinos/BUILD_DIR/Makefile.export.Trilinos_install to
ILU/epetra.

3. Edit the ILU/epetra/Makefile to set variables
a. MY_LIB_DIR to Trilinos/BUILD_DIR/lib 
b. RANLIB  to ranlib in your machine.

4. Do "make" or "make all" to compile HyperLU.

5. Edit ParaLU.xml to set the precondtioner (HyperLU/ILU/ILUT/ML currently) and
the matrix market file name. The matrix market file should be in the current 
directory.

5. mpirun -n np ./hyperlu_driver.exe will solve the linear system from the
matrix market file with a right hand side that is all ones.

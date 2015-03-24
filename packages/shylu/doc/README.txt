To Compile ShyLU
------------------

1. Compile Trilinos  with the following options
-DTrilinos_EXTRA_REPOSITORIES:STRING=preCopyrightTrilinos
-D Trilinos_ENABLE_ShyLU:BOOL=ON
-D Trilinos_ENABLE_TESTS:BOOL=ON (if you need the tests)

To test hyLU
------------------
1. You can use the two drivers in shylu/test,
AztecOO based driver (shylu_driver.exe) and a Belos based driver
(shylu_belos_driver.exe), to test your matrices. To use the different options
of ShyLU see step 2.

2. Edit ShyLU.xml to set the precondtioner (ShyLU)
and the matrix market file name.
    Please leave the parameter "Symmetry" as it is. I will end up doing a
    symmetric permutation of an unsymmetric matrix for now.

    The approximation method for the Schur complement is specified by
    "Schur Approximation Method". The current supported methods are
    "A22AndBlockDiagonals"  and "Threshold". "A22AndBlockDiagonals" uses the
    structure of A22 + a percentage of the block diagonals. Adjust the
    "Diagonal Factor" (0.0 to 1.0) for better convergence (at the cost of
    slower setup and more memory usage).
    If "Schur Approximation Method" is set to "Threshold" then ShyLU
    computes the entire Schur complement (not all at once) and drops entries
    based on "Relative Threshold". Setting "Relative Thresold" to 0.0 will
    compute the entire Schur complement.

    The inner solver is currently AztecOO with a domain decomposition 
    preconditioner. It is not exposed in the interface yet. I can do that if
    you need to change that.

3. mpirun -n np <driver_name> will solve the linear system from the
matrix market file with a right hand side that is all ones (or the one
you gave). This is probably
the way to run for now if you have the matrices in a matrix market file. 

To use ShyLU
------------------
1. To use ShyLU use Ifpack_ShyLU as your preconditioner. If
you use this interface you need to partition and redistribute the matrix
before you set the matrix in the preconditioner. Please see the driver for an
example. (This requirement is because of a known bug in AztecOO)

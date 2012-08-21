This directory contains some Matlab files to test 
a block partitioned Schur complement approach. 
SBBD.m computes the singly bordered block diagonal
form of a sparse matrix, but requires PaToH.
The LU*.m files supports factorization and solve of
a system in SBBD form. Currently the exact Schur complement
is computed. Future versions should compute approximate
Schur complements, as in ShyLU.

NOTE: The Matlab code is much less developed than
the ShyLU C++ code, and not intended for end users.
It may be useful for developers to prototype
new ideas or to generate sparsity plots.

#!/bin/csh
#
#  TestAmesos.csh, the Direct Sparse Solver Regresion Test, tests and times the 
#  direct sparse solvers supported by the AME interface.  At present, it 
#  tests two direct solvers:  Kundert and SuperLU MPI.  It also tests an 
#  indirect solver, AZTEC, for comparison.  
#
#  RETURNS 0 if the test succeeds and 1 if the test fails 
#
#  Each call to cxx_AME_mpi.exe performs one or more sovles on a particular 
#  matrix using a particular solver and prints a single line to SST.summary.
#  If the test worked, that line will contain OK.  Any line not containing 
#  the words OK, represents failure.
#  
#  Hence, the check for success is to make sure that all non-blank lines in 
#  SST.summary contain OK.
#
#  Each time that this script is run, it appends the result of any previous runs to 
#  AME.summary.  
#
#  More detailed loggin information can be found in SST.log
#
#  Planned enhancements:
#    0)  It should provde some output to reassure the user 
#    1)  Add time specifications - talk to Mike and Jim about how they do this
#    2)  Canwe put a time limit on the code?  
#    3)  Put the files into a subdirectory of $(TRILINOS_DATA_HOME) 
#    4)  Offer a -v option
#
#
#  A typical call to cxx_AME_mpi.exe is:
#COMMENT       mpirun -np 1 cxx_AME_mpi.exe SuperLUdist SuperLU.rua 0 1 1 0 1e-14 1e-14
#  where:
#     SuperLUdist SuperLU.rua - The solver to use and the matrix to solve
#     0 1 1 0                 - MatrixType, Special, NumSolves, Transpose
#     1e-14 1e-14             - max residual, max error 
#
#
#   MatType = 0 means serial (all stored on process 0) 
#   MatType = 1 means distributed (evenly) 
#   Special = 0 means use dgssvx (in SuperLU)
#   Special = 1 means use dgssv (in SuperLU)
#   NumSolves < 0 means use multiple right hand sides
#   NumSolves > 1 means use blocked right hand sides
#

cat >>AME.summary <SST.summary 
echo "COMMENT Start TestAmesos.csh, the Direct Sparse Solver Regresion Test" > SST.summary 
echo "COMMENT column 1 - machine name " >> SST.summary 
echo "COMMENT column 2 - solver name " >> SST.summary 
echo "COMMENT column 3 - timestamp" >> SST.summary 
echo "COMMENT column 4 - matrix file " >> SST.summary 
echo "COMMENT column 5 - Matrix type  " >> SST.summary 
echo "COMMENT column 6 - Special - only used for SuperLU serial " >> SST.summary 
echo "COMMENT column 7 - Number of processes " >> SST.summary 
echo "COMMENT column 8 - Number of right hand sides, nrhs, (-1 means multiple solves) " >> SST.summary 
echo "COMMENT column 9 - Tranpose (1 == solve A^t x = b)" >> SST.summary 
echo "COMMENT column 10 - Norm of the matrix " >> SST.summary 
echo "COMMENT column 11 - relative error - i.e. error/norm(X) " >> SST.summary 
echo "COMMENT column 12 - residual error - i.e. residual/norm(B) " >> SST.summary 
echo "COMMENT column 13 - total_time " >> SST.summary 
echo "COMMENT column 14 - MRHS only - Factor time " >> SST.summary 
echo "COMMENT column 15 - MRHS only - Time for one solve " >> SST.summary 
echo "COMMENT column 16 - MRHS only - Time for nrhs-1 solves " >> SST.summary 
echo "COMMENT column 17+ - summary " >> SST.summary 


#
#  Test one process, three processes and three processes transpose, tiny serial matrix, on SuperLUdist
#
yod -sz 1 cxx_AME_mpi.exe SuperLUdist SuperLU.rua 0 1 1 0 1e-14 1e-14 >>SST.stdout
yod -sz 3 cxx_AME_mpi.exe SuperLUdist SuperLU.rua 0 1 1 0 -2    -2  >>SST.stdout
yod -sz 3 cxx_AME_mpi.exe SuperLUdist SuperLU.rua 0 1 1 1 -2    -2  >>SST.stdout
#
#  Test one process, three processes and three processes transposes, tiny distributed matrix, on SuperLUdist
#
yod -sz 1 cxx_AME_mpi.exe SuperLUdist   fidapm05.rua 0 1 1 0 100000000 1 >>SST.stdout
yod -sz 3 cxx_AME_mpi.exe SuperLUdist   fidapm05.rua 1 1 1 0 100000000 1  >>SST.stdout
yod -sz 3 cxx_AME_mpi.exe SuperLUdist   fidapm05.rua 1 1 1 1 100000000 1  >>SST.stdout
#
#  Test some more small matrices
#
yod -sz 1 cxx_AME_mpi.exe SuperLUdist   ImpcolA.rua 0 1 1 0 1e-11 1e-12 >>SST.stdout
yod -sz 3 cxx_AME_mpi.exe SuperLUdist   ImpcolA.rua 0 1 1 0 1e-11 1e-12  >>SST.stdout
yod -sz 3 cxx_AME_mpi.exe SuperLUdist   ImpcolA.rua 0 1 1 1 1e-11 1e-12  >>SST.stdout
yod -sz 1 cxx_AME_mpi.exe SuperLUdist   ImpcolB.rua 0 1 1 0 1e-11 1e-13 >>SST.stdout
yod -sz 3 cxx_AME_mpi.exe SuperLUdist   ImpcolB.rua 0 1 1 0 1e-11 1e-14  >>SST.stdout
yod -sz 1 cxx_AME_mpi.exe SuperLUdist   ImpcolC.rua 0 1 1 0 1e-13 1e-13 >>SST.stdout
yod -sz 3 cxx_AME_mpi.exe SuperLUdist   ImpcolC.rua 0 1 1 0 1e-13 1e-13  >>SST.stdout
yod -sz 1 cxx_AME_mpi.exe SuperLUdist   ImpcolD.rua 0 1 1 0 1e-13 1e-13 >>SST.stdout
yod -sz 3 cxx_AME_mpi.exe SuperLUdist   ImpcolD.rua 0 1 1 0 1e-13 1e-13  >>SST.stdout
yod -sz 1 cxx_AME_mpi.exe SuperLUdist   ImpcolE.rua 0 1 1 0 1e-12 1e-11 >>SST.stdout
yod -sz 3 cxx_AME_mpi.exe SuperLUdist   ImpcolE.rua 0 1 1 0 1e-12 1e-11  >>SST.stdout
#
#  Test mid sized matrices on 1 and 4 processes, half of them starting out serial, 
#  half starting out distributed.  (OK, distributed has little meaning on one process, but ...)
#
yod -sz 1 cxx_AME_mpi.exe SuperLUdist   bcsstk24.rsa 0 1 1 0 1e-6  1e-1 >>SST.stdout
yod -sz 4 cxx_AME_mpi.exe SuperLUdist   bcsstk24.rsa 1 1 1 0 1e-6  1e-1 >>SST.stdout
yod -sz 1 cxx_AME_mpi.exe SuperLUdist   bcsstk18.rsa 1 1 1 0 1e-9 1e-4  >>SST.stdout
yod -sz 4 cxx_AME_mpi.exe SuperLUdist   bcsstk18.rsa 0 1 1 0 1e-9 1e-4  >>SST.stdout

#
#  Test some tranpose solves
#
yod -sz 4 cxx_AME_mpi.exe SuperLUdist   bcsstk18.rsa 0 1 1 1 1e-10 1e-4  >>SST.stdout
yod -sz 4 cxx_AME_mpi.exe SuperLUdist   bcsstk18.rsa 1 1 1 1 1e-10 1e-4  >>SST.stdout

#
#  Test blocked right hand sides
#
yod -sz 1 cxx_AME_mpi.exe SuperLUdist   ImpcolA.rua 0 1 2 0 1e-11 1e-12 >>SST.stdout
yod -sz 5 cxx_AME_mpi.exe SuperLUdist   ImpcolB.rua 0 1 4 0 1e-11 1e-14  >>SST.stdout
yod -sz 2 cxx_AME_mpi.exe SuperLUdist   ImpcolE.rua 0 1 6 0 1e-12 1e-11  >>SST.stdout
yod -sz 4 cxx_AME_mpi.exe SuperLUdist   bcsstk24.rsa 0 1 3 0 1e-6  1e-1 >>SST.stdout
yod -sz 1 cxx_AME_mpi.exe SuperLUdist   bcsstk18.rsa 0 1 5 0 1e-9 1e-4  >>SST.stdout
yod -sz 4 cxx_AME_mpi.exe SuperLUdist   bcsstk18.rsa 0 1 32 0 1e-9 1e-4  >>SST.stdout
#
#  Test multiple right hand sides
#
yod -sz 1 cxx_AME_mpi.exe SuperLUdist   ImpcolC.rua 0 1 -1 0 1e-13 1e-13 >>SST.stdout
yod -sz 5 cxx_AME_mpi.exe SuperLUdist   ImpcolD.rua 0 1 -2 0 1e-13 1e-13  >>SST.stdout
yod -sz 2 cxx_AME_mpi.exe SuperLUdist   ImpcolE.rua 0 1 -3 0 1e-12 1e-11  >>SST.stdout
yod -sz 4 cxx_AME_mpi.exe SuperLUdist   bcsstk24.rsa 0 1 -4 0 1e-6  1e-1 >>SST.stdout
yod -sz 1 cxx_AME_mpi.exe SuperLUdist   bcsstk18.rsa 0 1 -5 0 1e-9 1e-4  >>SST.stdout
yod -sz 4 cxx_AME_mpi.exe SuperLUdist   bcsstk18.rsa 0 1 -3 0 1e-9 1e-4  >>SST.stdout

#
#  Test some solves on Aztec
#
yod -sz 1 cxx_AME_mpi.exe AZTEC   SuperLU.rua  0 1 1 0 1e-14 1e-13 >>SST.stdout
yod -sz 1 cxx_AME_mpi.exe AZTEC   ImpcolA.rua  0 1 1 0 1e30 1e30 >>&SST.stdout
yod -sz 1 cxx_AME_mpi.exe AZTEC   bcsstk18.rsa 0 1 1 0 1e30 1e30  >>SST.stdout
yod -sz 1 cxx_AME_mpi.exe AZTEC   bcsstk24.rsa 1 1 1 0 1e30 1e30  >>SST.stdout


echo "\nCOMMENT End TestAmesos.csh" >> SST.summary 

#
#  Make sure that the tests ran 
#
set expected_lines = `grep mpirun TestAmesos.csh | grep -v COMMENT | wc`
set results = `grep OK SST.summary | wc`
if ($results[1] != $expected_lines[1] ) then
    echo 'I expected ' $expected_lines[1] ' correct test results, but only saw: ' $results[1] 
    grep -v OK SST.summary | grep -v COMMENT | grep " " && echo "Direct Sparse Solver Regression Test FAILED" 
    exit(1)
endif
#
#  Prints out success or failure and exit 
#
grep -v OK SST.summary | grep -v COMMENT | grep " " > /dev/null || echo "Direct Sparse Solver Regression Test passed on all" $expected_lines[1] " tests"
#
#  This should not generaly print anything as errors should have been caught in the if test above
#
grep -v OK SST.summary  | grep -v COMMENT | grep " " && echo "Direct Sparse Solver Regression Test FAILED" 
exit($status == 0)

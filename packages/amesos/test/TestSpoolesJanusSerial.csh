#!/bin/csh
#
#  TestSpoolesJanusSerial.csh, the Direct Sparse Solver Regresion Test, tests and times 
#  SpoolesOO() 
#
#  RETURNS 0 if the test succeeds and 1 if the test fails 
#
#  Each call to cxx_AME_mpi.exe performs one or more sovles on a particular 
#  matrix using SpoolesOO() and prints a single line to SST.summary.
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
#COMMENT       yod -sz 1 cxx_AME_mpi.exe SPOOLES SuperLU.rua 0 1 1 0 1e-14 1e-14
#  where:
#     SPOOLES SuperLU.rua - The solver to use and the matrix to solve
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
touch SST.summary
cat >>AME.summary <SST.summary 
echo "COMMENT Start TestSpoolesJanusSerial.csh, the Direct Sparse Solver Regresion Test" > SST.summary 
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
echo "COMMENT column 14 - multiple solves only - Factor time " >> SST.summary 
echo "COMMENT column 15 - multiple solves only - Time for one solve " >> SST.summary 
echo "COMMENT column 16 - multiple solves only - Time for nrhs-1 solves " >> SST.summary 
echo "COMMENT column 17+ - summary " >> SST.summary 


#
#  Test one process, three processes and three processes transpose, tiny serial matrix, on SPOOLES
#
yod -sz 1 cxx_AME_mpi.exe SPOOLES SuperLU.rua 0 1 1 0 1e-14 1e-14 >>SST.stdout
#
#  Test one process, three processes and three processes transposes, tiny distributed matrix, on SPOOLES
#
yod -sz 1 cxx_AME_mpi.exe SPOOLES   fidapm05.rua 0 1 1 0 1000 1e-12 >>SST.stdout
#
#  Test some more small matrices
#
yod -sz 1 cxx_AME_mpi.exe SPOOLES   ImpcolA.rua 0 1 1 0 1e-11 1e-12 >>SST.stdout
yod -sz 1 cxx_AME_mpi.exe SPOOLES   ImpcolB.rua 0 1 1 0 1e-10 1e-13 >>SST.stdout
yod -sz 1 cxx_AME_mpi.exe SPOOLES   ImpcolC.rua 0 1 1 0 1e-13 1e-13 >>SST.stdout
yod -sz 1 cxx_AME_mpi.exe SPOOLES   ImpcolD.rua 0 1 1 0 1e-12 1e-13 >>SST.stdout
yod -sz 1 cxx_AME_mpi.exe SPOOLES   ImpcolE.rua 0 1 1 0 1e-11 1e-11 >>SST.stdout
#
#  Test mid sized matrices on 1 and 4 processes, half of them starting out serial, 
#  half starting out distributed.  (OK, distributed has little meaning on one process, but ...)
#
#  NOTE:  Spooles does not handle distributed matrices at present
#
yod -sz 1 cxx_AME_mpi.exe SPOOLES   bcsstk24.rsa 0 1 1 0 1e-6  1e-1 >>SST.stdout
# Comment DISTRIBUTED MATRIX yod -sz 4 cxx_AME_mpi.exe SPOOLES   bcsstk24.rsa 1 1 1 0 1e-6  1e-1 >>SST.stdout
# Comment DISTRIBUTED MATRIX yod -sz 1 cxx_AME_mpi.exe SPOOLES   bcsstk18.rsa 1 1 1 0 1e-9 1e-4  >>SST.stdout

#
#  Test some tranpose solves
#
# Comment DISTRIBUTED MATRIX  yod -sz 4 cxx_AME_mpi.exe SPOOLES   ImpcolA.rua 1 1 1 1 1e-11 1e-12  >>SST.stdout


#
#  Test blocked right hand sides
#
# Comment Spooles does not handle BRHS yet  - mpirun  -np 1 cxx_AME_mpi.exe SPOOLES   ImpcolA.rua 0 1 2 0 1e-11 1e-12 >>SST.stdout
# Comment Spooles does not handle BRHS yet  - mpirun  -np 5 cxx_AME_mpi.exe SPOOLES   ImpcolB.rua 0 1 4 0 1e-11 1e-14  >>SST.stdout
# Comment Spooles does not handle BRHS yet  - mpirun  -np 2 cxx_AME_mpi.exe SPOOLES   ImpcolE.rua 0 1 6 0 1e-12 1e-11  >>SST.stdout
# Comment Spooles does not handle BRHS yet  - mpirun  -np 4 cxx_AME_mpi.exe SPOOLES   bcsstk24.rsa 0 1 3 0 1e-6  1e-1 >>SST.stdout
# Comment Spooles does not handle BRHS yet  - mpirun  -np 1 cxx_AME_mpi.exe SPOOLES   bcsstk18.rsa 0 1 5 0 1e-9 1e-4  >>SST.stdout
# Comment Spooles does not handle BRHS yet  - mpirun  -np 4 cxx_AME_mpi.exe SPOOLES   bcsstk18.rsa 0 1 32 0 1e-9 1e-4  >>SST.stdout
#
#  Test multiple right hand sides
#
# Comment Spooles does not handle MRHS yet - yod -sz 1 cxx_AME_mpi.exe SPOOLES   ImpcolC.rua 0 1 -1 0 1e-13 1e-13 >>SST.stdout
# Comment Spooles does not handle MRHS yet - yod -sz 5 cxx_AME_mpi.exe SPOOLES   ImpcolD.rua 0 1 -2 0 1e-13 1e-13  >>SST.stdout
# Comment Spooles does not handle MRHS yet - yod -sz 2 cxx_AME_mpi.exe SPOOLES   ImpcolE.rua 0 1 -3 0 1e-12 1e-11  >>SST.stdout
# Comment Spooles does not handle MRHS yet - yod -sz 4 cxx_AME_mpi.exe SPOOLES   bcsstk24.rsa 0 1 -4 0 1e-6  1e-1 >>SST.stdout
# Comment Spooles does not handle MRHS yet - yod -sz 1 cxx_AME_mpi.exe SPOOLES   bcsstk18.rsa 0 1 -5 0 1e-9 1e-4  >>SST.stdout
# Comment Spooles does not handle MRHS yet - yod -sz 4 cxx_AME_mpi.exe SPOOLES   bcsstk18.rsa 0 1 -3 0 1e-9 1e-4  >>SST.stdout

#
#  Test some triplet files
#  The .triU files are unsymmatric, the .triS files are symmetric, providing 
#  either the upper or lower triangular part.
#
yod -sz 1 cxx_AME_mpi.exe SPOOLES SuperLU.triU 0 1 1 0 1e-14 1e-14 >>SST.stdout

yod -sz 1 cxx_AME_mpi.exe SPOOLES K4989.triS 0 1 1 0 1e-10 1e-8 >>SST.stdout

yod -sz 1 cxx_AME_mpi.exe SPOOLES K5000.triS 0 1 1 0 1e-10 1e-8 >>SST.stdout


yod -sz 1 cxx_AME_mpi.exe SPOOLES Khead.triS 0 1 1 0 1e-13 1e-9 >>SST.stdout

#
#  Spooles completes, but returns a permuted answer for Khead.triS
#
yod -sz 4 cxx_AME_mpi.exe SPOOLES Khead.triS 0 1 1 0 -2    -2 >>SST.stdout


echo "\nCOMMENT End TestSpoolesJanusSerial.csh" >> SST.summary 

#
#  Make sure that the tests ran 
#
set expected_lines = `grep yod TestSpoolesJanusSerial.csh | grep -v COMMENT | wc`
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

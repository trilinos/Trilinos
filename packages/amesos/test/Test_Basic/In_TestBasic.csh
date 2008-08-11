#!/bin/csh
# ************************************************************************
# 
#                 Amesos: Direct Sparse Solver Package
#                 Copyright (2004) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************

#
#  BasicTest.exe, the Direct Sparse Solver Regresion Test, tests and times 
#
#  RETURNS 0 if the test succeeds and 1 if the test fails 
#
#  Each call to amesos_test.exe performs one or more sovles on a particular 
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
#  More detailed logging information can be found in SST.log
#
#  A typical call to amesos_test.exe is:
#COMMENT       $mpigo 1 ./amesos_test.exe KLU SuperLU.rua 0 1 1 0 1e-14 1e-13
#  where:
#     KLU SuperLU.rua - The solver to use and the matrix to solve
#     0 1 1 0                 - MatrixType, Special, NumSolves, Transpose
#     1e-14 1e-13             - max error, max residual
#
#
#   MatType = 0 means serial (all stored on process 0) 
#   MatType = 1 means distributed (evenly) 
#   Special = number of repeats (0 means run it once)
#   NumSolves < 0 means use multiple right hand sides
#   NumSolves > 1 means use blocked right hand sides
#
touch SST.summary
cat >>AME.summary <SST.summary 
echo "COMMENT Start BasicTest.exe, the Direct Sparse Solver Regresion Test" > SST.summary 
echo "COMMENT The values printed in columns 11 and 12 are relative." >> SST.summary
echo "COMMENT We test against absolute errors."   >> SST.summary
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

## Some machines use a command different than mpirun to run mpi jobs.  The
## test-harness.plx script sets the environment variable
## "TRILINOS_TEST_HARNESS_MPIGO_COMMAND".  We test for
## this value below.  If not set, we set it to a default value.

set mpigo = `printenv TRILINOS_TEST_HARNESS_MPIGO_COMMAND`

if ("$mpigo" == "") then
    set mpigo = "mpirun -np "
endif

if ("$1" == "") then
  echo "USAGE: TestBasic.csh [SOLVER]"
  echo "       SOLVER = KLU | UMFPACK | DSCPACK | SUPERLU | SUPERLUDIST | MUMPS"
  exit(1)
endif


# ====================== #
# non-symmetric matrices #
# ====================== #

if ("$1" != "DSCPACK" ) then

  # Test one process, three processes and three processes transposes, tiny distributed matrix
  $mpigo 1 ./amesos_test.exe $1   fidapm05.rua 0 1 1 0   1000000000000000 1e-1 >>SST.stdout
  $mpigo 2 ./amesos_test.exe $1   fidapm05.rua 1 1 4 1   1000000000000000 1e-1 >>SST.stdout
  $mpigo 4 ./amesos_test.exe $1   fidapm05.rua 1 1 -3 1  1000000000000000 1e-1 >>SST.stdout

  # Test some more small matrices
  $mpigo 1 ./amesos_test.exe $1   ImpcolA.rua 0 1 1 0  1e-9  1e-11  >>SST.stdout
  $mpigo 3 ./amesos_test.exe $1   ImpcolA.rua 0 0 1 0  1e-9  1e-11  >>SST.stdout
  $mpigo 3 ./amesos_test.exe $1   ImpcolA.rua 0 1 1 1  1e-9  1e-11  >>SST.stdout
  $mpigo 1 ./amesos_test.exe $1   ImpcolB.rua 0 1 1 0  1e-9  1e-12  >>SST.stdout
  $mpigo 3 ./amesos_test.exe $1   ImpcolB.rua 0 1 1 0  1e-9  1e-12  >>SST.stdout
  $mpigo 1 ./amesos_test.exe $1   ImpcolC.rua 0 1 1 0  1e-12 1e-13  >>SST.stdout
  $mpigo 3 ./amesos_test.exe $1   ImpcolC.rua 0 0 1 0  1e-12 1e-13  >>SST.stdout
  $mpigo 1 ./amesos_test.exe $1   ImpcolD.rua 0 1 1 0  1e-11 5e-13  >>SST.stdout
  $mpigo 3 ./amesos_test.exe $1   ImpcolD.rua 0 1 1 0  1e-11 5e-13  >>SST.stdout
  $mpigo 1 ./amesos_test.exe $1   ImpcolE.rua 0 1 1 0  1e-8  1e-9   >>SST.stdout
  $mpigo 3 ./amesos_test.exe $1   ImpcolE.rua 0 1 1 0  1e-8  1e-9   >>SST.stdout
  $mpigo 3 ./amesos_test.exe $1   ImpcolE.rua 0 1 3 1  1e-7  1e-9   >>SST.stdout
  $mpigo 3 ./amesos_test.exe $1   ImpcolE.rua 0 1 -3 1 1e-7  1e-9   >>SST.stdout
  $mpigo 3 ./amesos_test.exe $1   fidapm05.rua 1 1 1 1  1000000000000000 1e-1 >>SST.stdout

  #  Test some tranpose solves
  $mpigo 4 ./amesos_test.exe $1   ImpcolB.rua 0 1 1 1  1e-9 1e-12  >>SST.stdout
  $mpigo 4 ./amesos_test.exe $1   ImpcolA.rua 1 1 1 1  1e-9 1e-11  >>SST.stdout
  $mpigo 4 ./amesos_test.exe $1   ImpcolA.rua 1 1 3 1  1e-9 1e-11  >>SST.stdout
  $mpigo 4 ./amesos_test.exe $1   ImpcolA.rua 1 1 -2 1 1e-9 1e-11  >>SST.stdout

  #  Test blocked right hand sides
  $mpigo 1 ./amesos_test.exe $1   ImpcolA.rua 0 1 2 0 1e-9  1e-11  >>SST.stdout
  $mpigo 5 ./amesos_test.exe $1   ImpcolB.rua 0 1 4 0 1e-9  1e-12  >>SST.stdout
  $mpigo 2 ./amesos_test.exe $1   ImpcolE.rua 0 1 6 0 1e-7  1e-9   >>SST.stdout
  $mpigo 2 ./amesos_test.exe $1   ImpcolE.rua 0 1 6 1 1e-7  1e-9   >>SST.stdout

  #  Test multiple right hand sides
  $mpigo 1 ./amesos_test.exe $1   ImpcolC.rua 0 1 -1 0 1e-12 1e-13  >>SST.stdout
  $mpigo 5 ./amesos_test.exe $1   ImpcolD.rua 0 1 -2 0 1e-11 5e-13  >>SST.stdout
  $mpigo 5 ./amesos_test.exe $1   ImpcolD.rua 0 1 -2 1 1e-11 5e-13  >>SST.stdout
  $mpigo 2 ./amesos_test.exe $1   ImpcolE.rua 0 1 -3 0 1e-7  1e-9   >>SST.stdout

  #  Test blocked right hand sides with distributed matrix input
  $mpigo 1 ./amesos_test.exe $1   ImpcolA.rua 1 1 2 0 1e-9  1e-11  >>SST.stdout
  $mpigo 1 ./amesos_test.exe $1   ImpcolA.rua 1 1 2 1 1e-9  1e-11  >>SST.stdout
  $mpigo 5 ./amesos_test.exe $1   ImpcolB.rua 1 1 4 0 1e-9  1e-12  >>SST.stdout
  $mpigo 2 ./amesos_test.exe $1   ImpcolE.rua 1 1 6 0 1e-7  1e-9   >>SST.stdout

  #  Test multiple right hand sides with distributed matrix input
  $mpigo 1 ./amesos_test.exe $1   ImpcolC.rua 1 1 -2 0 1e-12 1e-13 >>SST.stdout
  $mpigo 5 ./amesos_test.exe $1   ImpcolD.rua 1 1 -3 1 1e-11 5e-13  >>SST.stdout
  $mpigo 5 ./amesos_test.exe $1   ImpcolD.rua 1 1 -3 0 1e-11 5e-13  >>SST.stdout
  $mpigo 2 ./amesos_test.exe $1   ImpcolE.rua 1 1 -1 0 1e-7  1e-10  >>SST.stdout

  #  A couple more matrices
  $mpigo 1 ./amesos_test.exe $1 662_bus_out.rsa 0 1 3 0 1e-8 1e-8 >>SST.stdout 

endif

# ====================== #
# now symmetric matrices #
# ====================== #

$mpigo 1 ./amesos_test.exe DSCPACK Diagonal.mtx 0 1 1 0 1e-15 1e-15 >>SST.stdout
$mpigo 1 ./amesos_test.exe DSCPACK bcsstk01.mtx 0 1 1 0 1e-12 1e-5 >>SST.stdout
$mpigo 1 ./amesos_test.exe DSCPACK bcsstk02.mtx 0 1 1 0 1e-12 1e-10 >>SST.stdout
$mpigo 1 ./amesos_test.exe DSCPACK bcsstk04.mtx 0 1 1 0 1e-11 1e-8 >>SST.stdout
$mpigo 1 ./amesos_test.exe DSCPACK nos1.mtx 0 1 1 0 1e-9 1e-5 >>SST.stdout
$mpigo 1 ./amesos_test.exe DSCPACK nos5.mtx 0 1 1 0 1e-12 1e-9 >>SST.stdout

$mpigo 2 ./amesos_test.exe DSCPACK Diagonal.mtx 0 1 1 0 1e-15 1e-15 >>SST.stdout
$mpigo 5 ./amesos_test.exe DSCPACK bcsstk01.mtx 0 1 1 0 1e-12 1e-5 >>SST.stdout
$mpigo 6 ./amesos_test.exe DSCPACK bcsstk02.mtx 0 1 1 0 1e-12 1e-10 >>SST.stdout
$mpigo 8 ./amesos_test.exe DSCPACK bcsstk04.mtx 0 1 1 0 1e-11 1e-8 >>SST.stdout
$mpigo 2 ./amesos_test.exe DSCPACK nos1.mtx 0 1 1 0 1e-9 1e-5 >>SST.stdout
$mpigo 2 ./amesos_test.exe DSCPACK nos5.mtx 0 1 1 0 1e-12 1e-9 >>SST.stdout

#  Test blocked right hand sides
$mpigo 4 ./amesos_test.exe DSCPACK bcsstk01.mtx 0 1 4 0 1e-12 1e-5 >>SST.stdout
$mpigo 6 ./amesos_test.exe DSCPACK bcsstk02.mtx 0 1 3 0 1e-12 1e-10 >>SST.stdout
$mpigo 3 ./amesos_test.exe DSCPACK bcsstk04.mtx 0 1 2 0 1e-11 1e-8 >>SST.stdout
$mpigo 4 ./amesos_test.exe DSCPACK nos1.mtx 0 1 2 0 1e-9 1e-5 >>SST.stdout
$mpigo 8 ./amesos_test.exe DSCPACK nos5.mtx 0 1 4 0 1e-12 1e-9 >>SST.stdout

#  Test multiple right hand sides
$mpigo 5 ./amesos_test.exe DSCPACK bcsstk01.mtx 0 1 -10 0 1e-12 1e-5 >>SST.stdout
$mpigo 6 ./amesos_test.exe DSCPACK bcsstk02.mtx 0 1 -3 0 1e-12 1e-10 >>SST.stdout
$mpigo 2 ./amesos_test.exe DSCPACK bcsstk04.mtx 0 1 -6 0 1e-11 1e-8 >>SST.stdout
$mpigo 2 ./amesos_test.exe DSCPACK nos1.mtx 0 1 -3 0 1e-9 1e-5 >>SST.stdout
$mpigo 2 ./amesos_test.exe DSCPACK nos5.mtx 0 1 -13 0 1e-12 1e-9 >>SST.stdout

echo "" >> SST.summary 
echo "COMMENT End BasicTest.exe" >> SST.summary 

# ============= #
# Check results #
# ============= #

set expected_lines = `grep -v COMMENT SST.summary | wc`
set results = `grep OK SST.summary | wc`
if ($results[1] != $expected_lines[1] - 1) then
    echo ' ' $expected_lines[1] - 1 ' tests ran, only: ' $results[1] ' succeeded '
    grep -v OK SST.summary | grep -v COMMENT | grep " " && echo "Direct Sparse Solver Regression Test FAILED" 
    exit(1)
endif

#  Prints out success or failure and exit 
grep -v OK SST.summary | grep -v COMMENT | grep " " > /dev/null || echo "Direct Sparse Solver Regression Test using $1 passed on all" $expected_lines[1] " tests that executed"

#  This should not generaly print anything as errors should have been caught in the if test above
grep -v OK SST.summary  | grep -v COMMENT | grep " " && echo "Direct Sparse Solver Regression Test using $1 FAILED" 
exit($status == 0)

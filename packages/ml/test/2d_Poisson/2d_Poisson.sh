#!/bin/sh
# ************************************************************************
# 
#               ML: A Multilevel Preconditioner Package
#                 Copyright (2002) Sandia Corporation
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
# Questions? Contact Jonathan Hu (jhu@sandia.gov)  or
#            Ray Tuminaro (rstumin@sandia.gov).
# 
# ************************************************************************
#
# This script is a wrapper for the executable "2D_poisson.exe".  It's assumed
# that "2D_poisson.exe" has been compiled with the macro ML_BENCHMARK defined.
#
# This script does the following:
#
# 1) In turn, each input deck "ml_inputfile[0-7]" is copied to "ml_inputfile".
#
# 2) "2D_poisson.exe" is run, and output is piped to "output[0-7]".
#
# 3) The results are compared to results saved in "baseline[0-7]".
#
# 4) A test is declared to have passed if either of the following are true:
#    a) the number of iterations is the same
#    b) wc reports no substantial differences in the results of diff applied
#       to the two files.
#    Otherwise, the test is said to have failed.
#
# Input
#    -v     Prints the file "output[0-7]" and "SUMMARY" to screen.
#
# ************************************************************************

#EXE="../../examples/ml_example2d.exe"
TESTNAME="2Dpoisson"
EXE="./ml_example2d.exe"
SUMMARY_FILE="SUMMARY"
CP="/bin/cp"
RM="/bin/rm"
VERBOSE="${1}tt"

${RM} -f output* difs* ml_inputfile

# execute any file that starts with ml_inputfile ...
for j in ml_inputfile*
do

#  j="ml_inputfile7"
  testnum=`echo $j | sed "s/ml_inputfile//"`
  ${CP} $j ml_inputfile

  # determine from input deck how many processes this job requires
  Nprocs=`grep -i "^Number of processors" ml_inputfile | sed s"/.*=//"`

  # default to one process if not specified in input deck
  temp=`echo xx${Nprocs}xx | sed "s/ *//g"`
  if test `expr ${temp}` = 'xxxx'
  then
      Nprocs=1
  fi

  # set the command to run (e.g., mpirun -np 4 ../../example/ml_example2d.exe)
  if test `expr ${Nprocs}` = '1'
  then
    EXE_MOD=$EXE
  else
    EXE_MOD="mpirun -np $Nprocs ${EXE}"
  fi

  # run
  ${EXE_MOD} >& output${testnum}

  if test `expr ${VERBOSE}` = '-vtt'
  then
    cat output${testnum}
  else
    diff -w baseline${testnum} output${testnum} > difs${testnum}
    # wc prints newline, word and byte counts
    echo "" >> ${SUMMARY_FILE}
    echo "**************************** regression test ${testnum} ******************************" >> ${SUMMARY_FILE}
    wc difs${testnum} | sed "s/difs${testnum}//" >> ${SUMMARY_FILE}
    # grep the dif file for significant changes:
    #   operator complexity, number of iterations, or solution entries
    cat difs${testnum} | grep "complexity" >> ${SUMMARY_FILE}
    check1="$?"
    cat difs${testnum} | grep "gid = " >> ${SUMMARY_FILE}
    check2="$?"
    grep "total iterations" difs${testnum} >> ${SUMMARY_FILE}
    check3="$?"


    # If any diffs occured, then the test fails.
    if (test "${check1}" == "0") || (test "${check2}" == "0") || (test "${check3}" == "0")
    then
      echo "********           ${TESTNAME}:    test #${testnum} FAILED"
    else
      echo "********           ${TESTNAME}:    test #${testnum} passed"
      ${RM} -f output${testnum}
    fi

  fi  #if test `expr ${VERBOSE}` = '-vtt'

  ${RM} -f difs${testnum}
  ${RM} -f ml_inputfile
  ${RM} -f output${testnum}
done

#!/bin/csh
set TESTNAME="2d_Poisson"
set EXE="./${TESTNAME}.exe"
set SUMMARY_FILE="SUMMARY"
set CP="/bin/cp"
set RM="/bin/rm"
set VERBOSE="${1}tt"
set EXITCODE=0

if ( ${VERBOSE} == '-vtt' ) then
  cat ${SUMMARY_FILE}
  ${RM} -f ${SUMMARY_FILE}
else
  # execute any file that starts with ml_inputfile ...
  foreach j ( ml_inputfile* )
  
    set testnum=`echo $j | sed "s/ml_inputfile//"`
    ${CP} $j ml_inputfile
  
    # determine from input deck how many processes this job requires
    set Nprocs=`grep -i "^Number of processors" ml_inputfile | sed s"/.*=//"`
  
    # default to one process if not specified in input deck
    set temp=`echo xx${Nprocs}xx | sed "s/ *//g"`
    if ( ${temp} == 'xxxx' ) then
        set Nprocs=1
    endif
  
    # set the command to run (e.g., mpirun -np 4 ../../example/ml_example2d.exe)
    if ( ${Nprocs} == '1' ) then
      set EXE_MOD=${EXE}
    else
      set EXE_MOD="mpirun -np $Nprocs ${EXE}"
    endif
  
    # run
    ${EXE_MOD} >& output${testnum}
  
    if ( ${VERBOSE} == '-vtt' ) then
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
      set check1="$?"
      cat difs${testnum} | grep "gid = " >> ${SUMMARY_FILE}
      set check2="$?"
      grep "total iterations" difs${testnum} >> ${SUMMARY_FILE}
      set check3="$?"
  
  
      # If any diffs occured, then the test fails.
      #if ( "${check3}" == "0" ) then
      if (("${check1}" == "0")||("${check2}" == "0")||("${check3}" == "0")) then
        echo "********           ${TESTNAME}:    test #${testnum} FAILED"  >> ${SUMMARY_FILE}
        set EXITCODE=1
      else
        echo "********           ${TESTNAME}:    test #${testnum} passed" >> ${SUMMARY_FILE}
        ${RM} -f output${testnum}
      endif
  
    ${RM} -f difs${testnum}
    ${RM} -f ml_inputfile
  end #foreach j ( ml_inputfile* )

  if ( ${EXITCODE} == 1 ) then
    foreach j ( output* )
      set testnum=`echo $j | sed "s/output//"`
      echo "" >> ${SUMMARY_FILE}
      echo "=============================================================================" >> ${SUMMARY_FILE}
      echo "==               output of failed regression test #${testnum}                 ==" >> ${SUMMARY_FILE}
      echo "=============================================================================" >> ${SUMMARY_FILE}
      cat ${j} >> ${SUMMARY_FILE}

      ${RM} -f ${j}
    end #foreach j ( ml_inputfile* )
  endif

endif  #if test `expr ${VERBOSE}` = '-vtt'

# 0 = all tests passed, 1 = at least one test failed
if ( ${EXITCODE} == 1 ) then
  echo " ***** Test ${TESTNAME} failed *****"
else
  ${RM} -f ${SUMMARY_FILE}
endif

#exit ${EXITCODE}

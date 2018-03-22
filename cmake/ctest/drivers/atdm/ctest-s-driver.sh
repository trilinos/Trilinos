#!/bin/bash -l

set +x

echo
echo "Start: ctest-s-driver.sh"
echo
echo "  ==> `date`"
echo

echo "Loading env and running ctest -S comamnd ..."

if [ "${WORKSPACE}" == ""  ] ; then
  echo "Error, must set WORKSPACE var before calling!"
  exit 1
fi

if [ "${JOB_NAME}" == ""  ] ; then
  echo "Error, must set JOB_NAME var before calling!"
  exit 1
fi

# Shiller/hansen settings for automated job (rest are set in ./environment.sh)
ulimit -c 0

export SUBDIR=SRC_AND_BUILD
if [ ! -e $SUBDIR ] ; then
  echo "Making $SUBDIR"
  mkdir $SUBDIR
fi

cd $SUBDIR/
echo "Current dir: $PWD"

source $WORKSPACE/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME
echo
module list

export Trilinos_REPOSITORY_LOCATION=https://github.com/trilinos/Trilinos.git

echo
echo "Running: ctest -V -S $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-driver.cmake ..."

ctest -V -S \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-driver.cmake
ATDM_TCD_CTEST_S_RETURN_CODE=$?

echo
echo "The 'ctest -S ctest-drivers.cmake' command returned code '$ATDM_TCD_CTEST_S_RETURN_CODE'"
echo
echo "End: ctest-s-driver.sh"
echo
echo "  ==> `date`"
echo
echo

# NOTE: The above output is important in order to determine if this script as
# run by the batch scheduler for the machine lets this script complete without
# killing it!

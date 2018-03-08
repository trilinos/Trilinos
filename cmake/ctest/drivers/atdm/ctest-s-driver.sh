#!/bin/bash -l

echo "Loading env and running ctest -S comamnd ..."

echo
echo "Time when the actual 'ctest -S ctest-driver.cmake' command was called:"
echo
echo "  ==> `date`"
echo

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

set -x

ctest -V -S \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-driver.cmake

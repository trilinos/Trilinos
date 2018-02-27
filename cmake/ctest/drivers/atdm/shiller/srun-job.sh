#!/bin/bash -l

if [ "${WORKSPACE}" == ""  ] ; then
  echo "Error, must set WORKSPACE var before calling"
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

source $WORKSPACE/Trilinos/cmake/std/atdm/shiller/environment.sh
echo
module list

set -x

ctest -V -S \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/shiller/ctest-driver.cmake

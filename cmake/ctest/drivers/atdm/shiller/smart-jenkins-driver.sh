#!/bin/bash -l
#
# This script is meant to be called directly from the Jenkins job.
#

if [ "${JOB_NAME}" == ""  ] ; then
  echo "Error, must set JOB_NAME var before calling!"
  exit 1
fi

if [ "${WORKSPACE}" == ""  ] ; then
  echo "Error, must set WORKSPACE var before calling"
  exit 1
fi

set -x

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/shiller/drivers/$JOB_NAME.sh

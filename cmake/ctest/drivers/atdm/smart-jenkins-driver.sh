#!/bin/bash -l
#
# This script is meant to be called directly from the Jenkins job.
#
# This script is not system-specific so it can be usesd on all systems.
#

echo "Running $JOB_NAME.sh script in smart Jenkins driver"

if [ "${JOB_NAME}" == ""  ] ; then
  echo "Error, must set JOB_NAME var before calling!"
  exit 1
fi

if [ "${WORKSPACE}" == ""  ] ; then
  echo "Error, must set WORKSPACE var before calling"
  exit 1
fi

source $WORKSPACE/Trilinos/cmake/std/atdm/utils/get_known_system_name.sh

echo
echo "Time when smart-jenkins-driver.sh was first called:"
echo
echo "  ==> `date`"
echo

set -x

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/$ATDM_CONFIG_KNOWN_SYSTEM_NAME/drivers/$JOB_NAME.sh

set +x

echo
echo "Time when smart-jenkins-driver.sh completed called:"
echo
echo "  ==> `date`"
echo
echo "DONE"
echo

#!/bin/bash -l
#
# This script is meant to be called directly from the Jenkins job.
#
# This script is not system-specific so it can be usesd on all systems.
#

echo "Start: smart-jenkins-drivers.sh"
echo
echo "  ==> `date`"
echo

if [ "${JOB_NAME}" == ""  ] ; then
  echo "Error, must set JOB_NAME var before calling!"
  exit 1
fi

if [ "${WORKSPACE}" == ""  ] ; then
  echo "Error, must set WORKSPACE var before calling"
  exit 1
fi

export ATDM_CONFIG_BUILD_NAME=$JOB_NAME
export ATDM_CONFIG_SCRIPT_DIR=$WORKSPACE/Trilinos/cmake/std/atdm
source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh

echo
echo "Running: $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/$ATDM_CONFIG_SYSTEM_NAME/drivers/$JOB_NAME.sh ..."

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/$ATDM_CONFIG_SYSTEM_NAME/drivers/$JOB_NAME.sh

echo
echo "End: smart-jenkins-drivers.sh"
echo
echo "  ==> `date`"
echo

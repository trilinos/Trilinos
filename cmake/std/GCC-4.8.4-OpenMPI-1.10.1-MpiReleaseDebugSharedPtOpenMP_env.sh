#
# Source this script in order to load env for
# GCC-4.8.4-MpiReleaseDebugSharedPtOpenMP build using:
#
#  source GCC-4.8.4-MpiReleaseDebugSharedPtOpenMP_env.sh
#

_ABS_FILE_PATH=`readlink -f "${BASH_SOURCE[0]}"`
if [ "$_ABS_FILE_PATH" != "" ] ; then
  _SCRIPT_DIR=`dirname $_ABS_FILE_PATH`
else
  echo "***"
  echo "*** Could not get the Trilinos directory!"
  echo "***"
  return
fi

source ${_SCRIPT_DIR}/../load_sems_dev_env.sh \
  sems-gcc/4.8.4 sems-openmpi/1.10.1 atdm-cmake/3.11.1

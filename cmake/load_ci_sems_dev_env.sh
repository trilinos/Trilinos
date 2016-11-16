# Load the standard Trilinos CI SEMS develoment environment
#
# This script is sourced as (takes no arugments):
#
#   $ source load_ci_sems_dev_env.sh
#
# All this script does is to source load_sems_dev_env.sh with the standard
# versions of the various modules used for CI testing.
#

# Get the base dir for the sourced script
called=$_
#[[ $called != $0 ]] && echo "Script is being sourced" || echo "Script is being run"
#echo "\$BASH_SOURCE ${BASH_SOURCE[@]}"
_SCRIPT_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
#echo "_SCRIPT_DIR = '$_SCRIPT_DIR'"

source $_SCRIPT_DIR/std/sems/get_default_modules.sh

if [ "$PLATFORM_TYPE" == "Linux" ] ; then
  DEFAULT_CI_SEMS_COMPILER=sems-gcc/4.7.2
elif [ "$PLATFORM_TYPE" == "Darwin" ] ; then
  DEFAULT_CI_SEMS_COMPILER=sems-gcc/5.3.0
else
  echo "ERROR, unknown platform type '$DEFAULT_CI_SEMS_COMPILER'!"
  exit 1
fi

source $_SCRIPT_DIR/load_sems_dev_env.sh \
  $DEFAULT_CI_SEMS_COMPILER sems-openmpi/1.6.5 sems-cmake/3.5.2

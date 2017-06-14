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

if [ "$1" != "" ] ; then
  echo "ERROR, the source script 'load_ci_sems_dev_env.sh' takes no arguments! (Remove '$1' ...)"
  return 1
fi

source $_SCRIPT_DIR/load_sems_dev_env.sh

# NOTE: Above, we will maintain the load_sems_dev_env.sh defaults so that they
# match the desired CI env.  That way, Trilinos users can accidentally source
# load_sems_dev_env.sh and they will get the stanard CI dev env.

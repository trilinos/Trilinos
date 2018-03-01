################################################################################
#
# Source to set up the env to do ATDM configuration of Trilinos.
#
################################################################################

# Assert this script is sourced, not run!
called=$_
if [ "$called" == "$0" ] ; then
  echo "This script '$0' is being called.  Instead, it must be sourced!"
  exit 1
fi

# Return the absoute directory of some relative directory path.
#
# This uses a temp shell to cd into the directory and then uses pwd to get the
# path.
function get_abs_dir_path() {
  [ -z "$1" ] && { pwd; return; }
  (cd -P -- "$1" && pwd)
}

# Get the base dir for the sourced script to find the base of Trilinos
_SCRIPT_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
#echo "_SCRIPT_DIR = '$_SCRIPT_DIR'"

#
# A) Parse the command-line arguments
#

# Make sure job-name is passed in as first (and only ) arguemnt
if [ "$1" == "" ] ; then
  echo "Error, the first argument must be the job name with keyword names!"
  return
fi

# Make sure there are no other command-line arguments set
if [ "$2" != "" ] ; then
  echo "Error, this source script only accepts a single comamnd-line argument!"
  return
fi

#
# B) Get the system name from the hostname
#

source $_SCRIPT_DIR/utils/get_known_system_name.sh

if [[ $ATDM_CONFIG_KNOWN_SYSTEM_NAME == "" ]] ; then
  echo "Error, could not determine known system, aborting env loading"
  return
fi

#
# C) Set JOB_NAME now that hostname has been asserted
#

echo "Setting export JOB_NAME=$1"
export JOB_NAME=$1

#
# D) Parse $JOB_NAME for consumption by the system-specific environoment.sh
# script
#

source $_SCRIPT_DIR/utils/set_build_options.sh

#
# E) Load the matching env
#

source $_SCRIPT_DIR/$ATDM_CONFIG_KNOWN_SYSTEM_NAME/environment.sh

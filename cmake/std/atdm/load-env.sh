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

THIS_HOSTNAME=`hostname`
#echo "Hostname = '$THIS_HOSTNAME'"

ATDM_HOSTNAME=

if [[ $THIS_HOSTNAME == "shiller"* ]] || [[ $THIS_HOSTNAME == "hansen"* ]] ; then
  ATDM_HOSTNAME=shiller
fi

# ToDo: Add more know hosts as you add them!

if [[ $ATDM_HOSTNAME == "" ]] ; then
  echo "Error, hostname = '$THIS_HOSTNAME' not recognized as a known ATDM system name!"
  return
else
  echo "Hostname '$THIS_HOSTNAME' matches known ATDM system '$ATDM_HOSTNAME'!"
fi

#
# C) Set JOB_NAME now that hostname has been asserted
#

echo "Setting JOB_NAME=$1"
export JOB_NAME=$1

#
# D) Load the matching env
#

source $_SCRIPT_DIR/$ATDM_HOSTNAME/environment.sh

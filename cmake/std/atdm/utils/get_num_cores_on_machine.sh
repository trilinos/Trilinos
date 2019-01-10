################################################################################
#
# Get the known system name (or error out)
#
# This script gets sourced from <system_name>/enviornment.sh scripts to
# determine how many cores are on a generic system.
#
################################################################################

unset ATDM_CONFIG_NUM_CORES_ON_MACHINE

# Assert this script is sourced, not run!
called=$_
if [ "$called" == "$0" ] ; then
  echo "This script '$0' is being called.  Instead, it must be sourced!"
  exit 1
fi
unset called

# Get num cores on rhel machines
export ATDM_CONFIG_NUM_CORES_ON_MACHINE=`grep processor /proc/cpuinfo | wc -l`

# Allow for override
if [[ "${ATDM_CONFIG_NUM_CORES_ON_MACHINE_OVERRIDE}" != "" ]] ; then
  export ATDM_CONFIG_NUM_CORES_ON_MACHINE=$ATDM_CONFIG_NUM_CORES_ON_MACHINE_OVERRIDE
fi

# ToDo: Add logic to find this out on other machines as well if needed!

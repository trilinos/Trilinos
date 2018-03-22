################################################################################
#
# Get the known system name (or error out)
#
################################################################################

# Assert this script is sourced, not run!
called=$_
if [ "$called" == "$0" ] ; then
  echo "This script '$0' is being called.  Instead, it must be sourced!"
  exit 1
fi

THIS_HOSTNAME=`hostname`
#echo "Hostname = '$THIS_HOSTNAME'"

ATDM_HOSTNAME=
ATDM_SYSTEM_NAME=

if [[ $THIS_HOSTNAME == "hansen"* ]] ; then
  ATDM_HOSTNAME=hansen
  ATDM_SYSTEM_NAME=shiller
elif [[ $THIS_HOSTNAME == "shiller"* ]] ; then
  ATDM_HOSTNAME=shiller
  ATDM_SYSTEM_NAME=shiller
elif [[ $THIS_HOSTNAME == "white"* ]] ; then
  ATDM_HOSTNAME=white
  ATDM_SYSTEM_NAME=ride
elif [[ $THIS_HOSTNAME == "ride"* ]] ; then
  ATDM_HOSTNAME=ride
  ATDM_SYSTEM_NAME=ride
fi

# ToDo: Add more know hosts as you add them!

if [[ $ATDM_SYSTEM_NAME == "" ]] ; then
  echo "Error, hostname = '$THIS_HOSTNAME' not recognized as a known ATDM system name!"
  return
else
  echo "Hostname '$THIS_HOSTNAME' matches known ATDM host '$ATDM_HOSTNAME' and system '$ATDM_SYSTEM_NAME'"
fi

export ATDM_CONFIG_KNOWN_HOSTNAME=$ATDM_HOSTNAME
export ATDM_CONFIG_KNOWN_SYSTEM_NAME=$ATDM_SYSTEM_NAME

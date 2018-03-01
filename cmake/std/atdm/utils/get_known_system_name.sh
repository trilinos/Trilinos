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

if [[ $THIS_HOSTNAME == "shiller"* ]] || [[ $THIS_HOSTNAME == "hansen"* ]] ; then
  ATDM_HOSTNAME=shiller
fi

# ToDo: Add more know hosts as you add them!

if [[ $ATDM_HOSTNAME == "" ]] ; then
  echo "Error, hostname = '$THIS_HOSTNAME' not recognized as a known ATDM system name!"
  return
else
  echo "Hostname '$THIS_HOSTNAME' matches known ATDM system '$ATDM_HOSTNAME'"
fi

export ATDM_CONFIG_KNOWN_SYSTEM_NAME=$ATDM_HOSTNAME

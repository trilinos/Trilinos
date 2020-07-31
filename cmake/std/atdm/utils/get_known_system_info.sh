################################################################################
#
# Get the known system name (or error out)
#
################################################################################

# Clean out vars in case this crashes before finishing

unset ATDM_CONFIG_REAL_HOSTNAME
unset ATDM_CONFIG_CDASH_HOSTNAME
unset ATDM_CONFIG_SYSTEM_NAME
unset ATDM_CONFIG_SYSTEM_DIR

called=$_
if [ "$called" == "$0" ] ; then
  echo "This script '$0' is being called.  Instead, it must be sourced!"
  exit 1
fi
unset called

if [ -z "$ATDM_CONFIG_BUILD_NAME" ] ; then
  echo "Error, must set ATDM_CONFIG_BUILD_NAME in env!"
  return
fi

if [ -z "$ATDM_CONFIG_SCRIPT_DIR" ] ; then
  echo "Error, must set ATDM_CONFIG_SCRIPT_DIR in env!"
  return
fi

source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info_utils.sh

if [[ "${ATDM_CONFIG_GET_KNOW_SYSTEM_INFO_REAL_HOSTNAME_OVERRIDE_FOR_UNIT_TESTING}" ]] ; then
  if [[ -z $ATDM_CONFIG_DISABLE_WARNINGS ]]; then
    echo
    echo "***"
    echo "*** WARNING: realHostname=$realHostname overriden to value of"
    echo "*** ATDM_CONFIG_GET_KNOW_SYSTEM_INFO_REAL_HOSTNAME_OVERRIDE_FOR_UNIT_TESTING='${ATDM_CONFIG_GET_KNOW_SYSTEM_INFO_REAL_HOSTNAME_OVERRIDE_FOR_UNIT_TESTING}'"
    echo "*** in <trilinos-dir>/cmake/std/atdm/utils/get_known_system_info.sh."
    echo "*** This variable should only be set for unit testing purposes!"
    echo "***"
    echo
  fi
  realHostname=${ATDM_CONFIG_GET_KNOW_SYSTEM_INFO_REAL_HOSTNAME_OVERRIDE_FOR_UNIT_TESTING}
else
  realHostname=`hostname`
  ATDM_CONFIG_SEMS_GET_PLATFORM=/projects/sems/modulefiles/utils/get-platform
fi
#echo "Hostname = '$realHostname'"

#
# List out all of the known system envs
#
# The order these are listed in this array matters only if multiple known
# system name keywords are listed in the build name string.  For example, if
# both 'ride' and 'cts1' are listed in the build name, then 'ride' will be the
# one recognized and 'cts1' will be ignored (because 'ride' is listed above
# 'cts1').
#
# However, it is important that "all" of the known systems be listed in this
# array in order for each system name to recognized in the build name string.
#

ATDM_KNOWN_SYSTEM_NAMES_LIST=(
  shiller
  ride
  ats1
  mutrino   # Deprecated, to be repalced by 'ats1'
  waterman
  ats2
  van1-tx2
  cts1empire
  cts1
  tlcc2
  sems-rhel7
  sems-rhel6
  cee-rhel6  # Used for CEE RHEL7 machines as well!
  spack-rhel
  )

#
# A) Look for matches of known system names that appear in the buildname
#

knownSystemNameInBuildName=`get_knownSystemNameInBuildName`
#echo "knownSystemNameInBuildName = '${knownSystemNameInBuildName}'"

# System name and hostname matching
systemNameTypeMatchedList=()  # In order of match preference
unset systemNameTypeMatchedListHostNames
declare -A systemNameTypeMatchedListHostNames

#
# B) See if the current system matches a known hostname
#
# This will be just a single match (if any)
#

hostnameMatch=
hostnameMatchSystemName=

# Specifically named test-bed machines
if [[ $realHostname == "hansen"* ]] ; then
  hostnameMatch=hansen
  hostnameMatchSystemName=shiller
elif [[ $realHostname == "shiller"* ]] ; then
  hostnameMatch=shiller
  hostnameMatchSystemName=shiller
elif [[ $realHostname == "white"* ]] ; then
  hostnameMatch=white
  hostnameMatchSystemName=ride
elif [[ $realHostname == "ride"* ]] ; then
  hostnameMatch=ride
  hostnameMatchSystemName=ride
elif [[ $realHostname == "waterman"* ]] ; then
  hostnameMatch=waterman
  hostnameMatchSystemName=waterman
elif [[ $realHostname == "vortex"* ]] ; then
  hostnameMatch=vortex
  hostnameMatchSystemName=ats2
# End specifically named systems
fi

#echo "hostnameMatch ='${hostnameMatch}'"

if [[ "${hostnameMatch}" != "" ]] ; then
  # A matching system by hostname becomes the first preferred match (but not
  # the only possible match)
  systemNameTypeMatchedList+=(${hostnameMatchSystemName})
  systemNameTypeMatchedListHostNames[${hostnameMatchSystemName}]=${hostnameMatch}
fi

#
# C) Look for known system types that matches this machine
#
# A given machine may match more than one known system type so this can be an
# array of matches.  Also, the list of known systems will be in the prefered
# match order so, if no other match criteria is in play, then the first
# matching system type will be selected.
#

# ATS-1 systems
if [[ $realHostname == "mutrino"* || $HOST == "mutrino"* ]] ; then
  systemNameTypeMatchedList+=(ats1)
  systemNameTypeMatchedListHostNames[ats1]=mutrino
  systemNameTypeMatchedList+=(mutrino)
  systemNameTypeMatchedListHostNames[mutrino]=mutrino
fi

# ASTRA/Van1-Tx2 systems
if [[ $SNLSYSTEM == "astra"* ]] ; then
  systemNameTypeMatchedList+=(van1-tx2)
  systemNameTypeMatchedListHostNames[van1-tx2]=$SNLCLUSTER
fi

# CTS1 systems
if [[ $SNLSYSTEM == "cts1" ]] ; then
  # Make cts1empire the default environment
  systemNameTypeMatchedList+=(cts1empire)
  systemNameTypeMatchedListHostNames[cts1empire]=$SNLCLUSTER
  # Add cts1 to the list for the D.1 branch, below
  systemNameTypeMatchedList+=(cts1)
  systemNameTypeMatchedListHostNames[cts1]=$SNLCLUSTER
fi

# TLCC2 systems
if [[ $SNLSYSTEM == "tlcc2"* ]] ; then
  systemNameTypeMatchedList+=(tlcc2)
  systemNameTypeMatchedListHostNames[tlcc2]=$SNLCLUSTER
fi

# SEMS RHEL6 and RHEL7 systems
if [[ "${SEMS_PLATFORM}" == "rhel6-x86_64" ]] ; then
  systemNameTypeMatchedList+=(sems-rhel6)
  systemNameTypeMatchedListHostNames[sems-rhel6]=sems-rhel6
elif [[ "${SEMS_PLATFORM}" == "rhel7-x86_64" ]] ; then
  systemNameTypeMatchedList+=(sems-rhel7)
  systemNameTypeMatchedListHostNames[sems-rhel7]=sems-rhel7
elif [[ "${SNLSYSTEM}" == "astra" || \
        "${SNLSYSTEM}" == "vortex" ]] ; then
  echo "Don't call get-platform on 'astra' systems" > /dev/null
  # Above logic avoids an 'ERROR: Unrecognized cluster <name>' on these systems
elif [[ -f $ATDM_CONFIG_SEMS_GET_PLATFORM ]] ; then
  ATDM_SYSTEM_NAME=`source $ATDM_CONFIG_SEMS_GET_PLATFORM`
  if [[ $ATDM_SYSTEM_NAME == "rhel6-x86_64" ]] ; then
    systemNameTypeMatchedList+=(sems-rhel6)
    systemNameTypeMatchedListHostNames[sems-rhel6]=sems-rhel6
  elif [[ $ATDM_SYSTEM_NAME == "rhel7-x86_64" ]] ; then
    systemNameTypeMatchedList+=(sems-rhel7)
    systemNameTypeMatchedListHostNames[sems-rhel7]=sems-rhel7
  fi
fi

# CEE RHEL6 (and RHEL7) systems
if [[ "${SNLSYSTEM}" == "cee" ]] ; then
  if [[ "${SNLCLUSTER}" == "linux_rh6" ]] || [[ "${SNLCLUSTER}" == "linux_rh7" ]] ; then
    systemNameTypeMatchedList+=(cee-rhel6)
    systemNameTypeMatchedListHostNames[cee-rhel6]=cee-rhel6
  fi
fi

# If the user puts 'spack-rhel' in the build name, assume that the modules are
# there (since one can build this stack on almost any machine).
if [[ "${knownSystemNameInBuildName}" == "spack-rhel" ]] ; then
  systemNameTypeMatchedList+=(spack-rhel)
  systemNameTypeMatchedListHostNames[spack-rhel]=spack-rhel
fi
# NOTE: Above, I was using:
#
#  spackCMakeModule=`module avail 2>&1 | grep spack-cmake | head -1`
#
# to determine if a spack-rhel env was setup but that check actually takes a
# long time (5 sec) on some systems so I went with the easier logic above.

#echo "systemNameTypeMatchedList = '(${systemNameTypeMatchedList[@]})'"

#################################################################################
### NOTE: No modifications below this line should be needed when adding a new
### system base on hostname or system type!
#################################################################################

#
# D) Select a known system given the above info
#

ATDM_HOSTNAME=
ATDM_SYSTEM_NAME=

# D.1) First, go with the system name in the build name if one was recognised
if [[ "${ATDM_SYSTEM_NAME}" == "" ]] && [[ "${knownSystemNameInBuildName}" != "" ]] ; then
  ATDM_SYSTEM_NAME=${knownSystemNameInBuildName}
  ATDM_HOSTNAME=${systemNameTypeMatchedListHostNames[${ATDM_SYSTEM_NAME}]}
  assert_selected_system_matches_known_system_type_matches || return
fi

# D.2) Last, go with the first matching system name on this machine
if [[ "${ATDM_SYSTEM_NAME}" == "" ]] && [[ "${systemNameTypeMatchedList}" != "" ]] ; then
  ATDM_SYSTEM_NAME=${systemNameTypeMatchedList[0]}  # First matching system type is preferred!
  ATDM_HOSTNAME=${systemNameTypeMatchedListHostNames[${ATDM_SYSTEM_NAME}]}
fi

#echo "ATDM_HOSTNAME = '${ATDM_HOSTNAME}'"
#echo "ATDM_SYSTEM_NAME = '${ATDM_SYSTEM_NAME}'"

#
# E) We have selected a known system set the env vars for that!
#

if [[ $ATDM_SYSTEM_NAME != "" ]] ; then
  echo "Hostname '$realHostname' matches known ATDM host '$realHostname' and system '$ATDM_SYSTEM_NAME'"
  export ATDM_CONFIG_REAL_HOSTNAME=$realHostname
  export ATDM_CONFIG_CDASH_HOSTNAME=$ATDM_HOSTNAME
  export ATDM_CONFIG_SYSTEM_NAME=$ATDM_SYSTEM_NAME
  export ATDM_CONFIG_SYSTEM_DIR=${ATDM_CONFIG_SCRIPT_DIR}/${ATDM_CONFIG_SYSTEM_NAME}
fi

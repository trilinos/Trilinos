################################################################################
#
# Get the known system name (or error out)
#
################################################################################

# Clean out vars in case this crashes before finishing

unset ATDM_CONFIG_REAL_HOSTNAME
unset ATDM_CONFIG_KNOWN_HOSTNAME
unset ATDM_CONFIG_KNOWN_SYSTEM_NAME

# Assert this script is sourced, not run!
called=$_
if [ "$called" == "$0" ] ; then
  echo "This script '$0' is being called.  Instead, it must be sourced!"
  exit 1
fi
unset called

# Assert that ATDM_CONFIG_BUILD_NAME is set!
if [ -z "$ATDM_CONFIG_BUILD_NAME" ] ; then
  echo "Error, must set ATDM_CONFIG_BUILD_NAME in env!"
  return
fi

export ATDM_CONFIG_REAL_HOSTNAME=`hostname`
#echo "Hostname = '$ATDM_CONFIG_REAL_HOSTNAME'"

ATDM_HOSTNAME=
ATDM_SYSTEM_NAME=

ATDM_IS_CEE_RHEL6_MACHINE=

if [[ $ATDM_CONFIG_REAL_HOSTNAME == "hansen"* ]] ; then
  ATDM_HOSTNAME=hansen
  ATDM_SYSTEM_NAME=shiller
elif [[ $ATDM_CONFIG_REAL_HOSTNAME == "shiller"* ]] ; then
  ATDM_HOSTNAME=shiller
  ATDM_SYSTEM_NAME=shiller
elif [[ $ATDM_CONFIG_REAL_HOSTNAME == "white"* ]] ; then
  ATDM_HOSTNAME=white
  ATDM_SYSTEM_NAME=ride
elif [[ $ATDM_CONFIG_REAL_HOSTNAME == "ride"* ]] ; then
  ATDM_HOSTNAME=ride
  ATDM_SYSTEM_NAME=ride
elif [[ $ATDM_CONFIG_REAL_HOSTNAME == "chama"* ]] ; then
  ATDM_HOSTNAME=chama
  ATDM_SYSTEM_NAME=chama
elif [[ $ATDM_CONFIG_REAL_HOSTNAME == "serrano"* ]] \
  || [[ $ATDM_CONFIG_REAL_HOSTNAME =~ ser[0-9]+ ]] ; then
  ATDM_HOSTNAME=serrano
  ATDM_SYSTEM_NAME=serrano
elif [[ $ATDM_CONFIG_REAL_HOSTNAME == "eclipse"* ]] \
  || [[ $ATDM_CONFIG_REAL_HOSTNAME =~ ec[0-9]+ ]] ; then
  ATDM_HOSTNAME=eclipse
  ATDM_SYSTEM_NAME=serrano
elif [[ $ATDM_CONFIG_REAL_HOSTNAME == "ghost"* ]] \
  || [[ $ATDM_CONFIG_REAL_HOSTNAME =~ gho[0-9]+ ]] ; then
  ATDM_HOSTNAME=ghost
  ATDM_SYSTEM_NAME=serrano
elif [[ $ATDM_CONFIG_REAL_HOSTNAME == "mutrino"* ]] ; then
  ATDM_HOSTNAME=mutrino
  ATDM_SYSTEM_NAME=mutrino
elif [[ $ATDM_CONFIG_REAL_HOSTNAME == "waterman"* ]] ; then
  ATDM_HOSTNAME=waterman
  ATDM_SYSTEM_NAME=waterman
elif [[ -f /projects/sems/modulefiles/utils/get-platform ]] ; then
  # This machine has the SEMS modules!
  ATDM_SYSTEM_NAME=`source /projects/sems/modulefiles/utils/get-platform`
  if [[ $ATDM_SYSTEM_NAME == "rhel7-x86_64" ]] ; then
    # This is a RHEL7 platform that has the SEMS modules.
    ATDM_HOSTNAME=sems-rhel7
    ATDM_SYSTEM_NAME=sems-rhel7
  elif [[ $ATDM_SYSTEM_NAME == "rhel6-x86_64" ]] ; then
    # This is a RHEL6 platform that has the SEMS modules.  But is this also a
    # CEE LAN mahcine?
    if [[ -f /projects/sparc/modules/cee-rhel6/sparc/master ]] ; then
      ATDM_IS_CEE_RHEL6_MACHINE=1
    fi
    # Now select the env based on the above logic
    if [[ $ATDM_CONFIG_BUILD_NAME == *"sems-rhel6"* ]] ; then
      ATDM_HOSTNAME=sems-rhel6
      ATDM_SYSTEM_NAME=sems-rhel6
    elif [[ $ATDM_CONFIG_BUILD_NAME == *"cee-rhel6"* ]] ; then
      if [[ $ATDM_IS_CEE_RHEL6_MACHINE == "1" ]] ; then
        # This is a CEE RHEL6 machine and 'cee-rhel6' was given in build name,
        # so select the system name 'sem-rhel6'
        ATDM_SYSTEM_NAME=cee-rhel6
        ATDM_HOSTNAME=cee-rhel6
      else
        echo
        echo "***"
        echo "*** Error, hostname='$ATDM_CONFIG_REAL_HOSTNAME' is a 'sems-rhel6' machine but is"
        echo "*** is *not* a 'cee-rhel6' machine but 'cee-rhel6' was given in the"
        echo "*** build name string '$ATDM_CONFIG_BUILD_NAME'!  Please remove 'cee-rhel6'"
	echo "*** from the build name or provide 'sems-rhel6' and then the 'sems-rhel6'"
	echo "*** env will be used."
        echo "***"
        return
      fi
    else
      # 'cee-rhel6' nor 'sems-rhel6' was set in the build name, so the default
      # is to sue the 'sems-rhel6' env!
      ATDM_HOSTNAME=sems-rhel6
      ATDM_SYSTEM_NAME=sems-rhel6
    fi
  else
    echo
    echo "***"
    echo "*** Error, hostname='$ATDM_CONFIG_REAL_HOSTNAME' has the SEMS env"
    echo "*** mounted but the SEMS system '${ATDM_SYSTEM_NAME}' is not yet supported!"
    echo "***"
    return
  fi
fi

if [[ $ATDM_SYSTEM_NAME == "" ]] ; then
  echo
  echo "***"
  echo "*** Error, hostname = '$ATDM_CONFIG_REAL_HOSTNAME' not recognized as a known ATDM system name!"
  echo "***"
  return
else
  echo "Hostname '$ATDM_CONFIG_REAL_HOSTNAME' matches known ATDM host '$ATDM_HOSTNAME' and system '$ATDM_SYSTEM_NAME'"
fi

export ATDM_CONFIG_KNOWN_HOSTNAME=$ATDM_HOSTNAME
export ATDM_CONFIG_KNOWN_SYSTEM_NAME=$ATDM_SYSTEM_NAME

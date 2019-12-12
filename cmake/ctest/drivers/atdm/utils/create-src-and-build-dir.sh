#
# Helper script that creates the SRC_AND_BUILD subdir and moves in
#

if [ "${WORKSPACE}" == ""  ] ; then
  echo "Error, must set WORKSPACE var before calling!"
  exit 1
fi

if [ "${JOB_NAME}" == ""  ] ; then
  echo "Error, must set JOB_NAME var before calling!"
  exit 1
fi

# Don't allow core file dumps (which take a huge amount of memory for lots of
# failing tests)
ulimit -c 0

# Set up the SRC_AND_BUILD subdir

if [[ "${ATDM_CONFIG_WORKSPACE_BASE}" != "" ]] ; then

  ATDM_CONFIG_WORKSPACE="${ATDM_CONFIG_WORKSPACE_BASE}/${ATDM_CONFIG_SYSTEM_NAME}/${JOB_NAME}"

  echo
  echo "Using alternate workspace dir '${ATDM_CONFIG_WORKSPACE}'" 

  if [[ ! -e "${ATDM_CONFIG_WORKSPACE}" ]] ; then
    echo
    echo "Creating non-existent workspace dir '${ATDM_CONFIG_WORKSPACE}'"
    mkdir -p "${ATDM_CONFIG_WORKSPACE}"
  fi

  ATDM_CONFIG_WORKSPACE_SRC_AND_BUILD="${ATDM_CONFIG_WORKSPACE}/SRC_AND_BUILD"

  # First copy an existing SRC_AND_BUILD directory to transition to the new
  # workspace base dir
  if [[ -e SRC_AND_BUILD ]] && [[ ! -e "${ATDM_CONFIG_WORKSPACE_SRC_AND_BUILD}" ]] ; then
    echo
    echo "Local 'SRC_AND_BUILD' does exist but '${ATDM_CONFIG_WORKSPACE_SRC_AND_BUILD}' does not so move it there (to preserve cloned Trilinos repo)"
    mv SRC_AND_BUILD "${ATDM_CONFIG_WORKSPACE_SRC_AND_BUILD}"
  fi

  # Create target SRC_AND_BUILD dir if does not exist
  if [[ ! -e "${ATDM_CONFIG_WORKSPACE_SRC_AND_BUILD}" ]] ; then
    echo
    echo "Create new workspace dir '${ATDM_CONFIG_WORKSPACE_SRC_AND_BUILD}'"
    mkdir "${ATDM_CONFIG_WORKSPACE_SRC_AND_BUILD}"
  fi

  # Set up local SRC_AND_BUILD symlink to alternate workspace
  if [[ -e SRC_AND_BUILD ]] ; then
    if [[ `atdm_config_get_abs_dir_path SRC_AND_BUILD` == "${ATDM_CONFIG_WORKSPACE_SRC_AND_BUILD}" ]] ; then
      echo
      echo "Local SRC_AND_BUILD already symlinked to correct alternate workspace '${ATDM_CONFIG_WORKSPACE_SRC_AND_BUILD}'!"
    else
      echo
      echo "Local SRC_AND_BUILD not correctly symlinked to correct alternate workspace so remove it!"
      rm -rf SRC_AND_BUILD
    fi
  fi

  # Create symlink for SRC_AND_BUILD if not already set up
  if [[ ! -e SRC_AND_BUILD ]] ; then
    echo
    echo "Creating local symlink to alternate workspace '${ATDM_CONFIG_WORKSPACE_SRC_AND_BUILD}'!"
    ln -s "${ATDM_CONFIG_WORKSPACE_SRC_AND_BUILD}" SRC_AND_BUILD
  fi

  # Change the workspace to actually be the alternate directory (critical for the install)
  cd "${ATDM_CONFIG_WORKSPACE}"
  ln -s "${WORKSPACE}/Trilinos" .  # Must have base Trilinos dir like Jenkins cloned it
  cd -
  export WORKSPACE="${ATDM_CONFIG_WORKSPACE}"

else

  # We are using local workspace directly under this current directory.

  if [ ! -e SRC_AND_BUILD ] ; then
    echo "Making local SRC_AND_BUILD"
    mkdir SRC_AND_BUILD
  fi

fi

# Move into the correct workspace subdir

cd "${WORKSPACE}/SRC_AND_BUILD/"
echo
echo "Current dir: $PWD"


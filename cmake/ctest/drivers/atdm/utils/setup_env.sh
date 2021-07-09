set +x

#
# A) Load the env
#

if [[ "${Trilinos_ENABLE_BUILD_STATS}" == "" ]] ; then
  export Trilinos_ENABLE_BUILD_STATS=OFF
fi
echo "Trilinos_ENABLE_BUILD_STATS='${Trilinos_ENABLE_BUILD_STATS}'"

source ${WORKSPACE}/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME
echo
module list
echo
echo "cmake in path:"
which cmake
echo
if [[ "$ATDM_CONFIG_USE_NINJA" == "ON" ]]; then
  echo "ninja in path:"
  which ninja
fi
echo
echo "ATDM config env vars:"
set | grep ATDM_CONFIG_
echo
echo "OpenMP env vars:"
set | grep ^OMP_
echo
echo "PATH=$PATH"

#
# B) Set the defaut git repo for Trilinos
#

if [ "${Trilinos_REPOSITORY_LOCATION}" == "" ] ; then
  export Trilinos_REPOSITORY_LOCATION=https://github.com/trilinos/Trilinos.git
fi

unset http_proxy
unset HTTP_PROXY
unset no_proxy
unset NO_PROXY
# NOTE: Above we have to unset http_proxy and no_proxy to allow the
# cdash submits to go  through. Note that we can't  unset https_proxy
# which is needed for the git operations with  https://github.com.

#
# C) Setup install-releated stuff
#

source ${ATDM_CONFIG_SCRIPT_DIR}/atdm_devops_install_defaults.sh

echo

# Set up default ATDM_CONFIG_USE_XXX_DEFAULT vars from
# ATDM_CONFIG_USE_JENKINS_INSTALL_DEFAULTS

if  [[ "${ATDM_CONFIG_USE_WORKSPACE_BASE_DEFAULT}" == "" ]] \
  && [[ "${ATDM_CONFIG_USE_JENKINS_INSTALL_DEFAULTS}" == "1" ]] ; then
  export ATDM_CONFIG_USE_WORKSPACE_BASE_DEFAULT=1
fi
echo "ATDM_CONFIG_USE_WORKSPACE_BASE_DEFAULT=${ATDM_CONFIG_USE_WORKSPACE_BASE_DEFAULT}"

if  [[ "${ATDM_CONFIG_USE_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE_DEFAULT}" == "" ]] \
  && [[ "${ATDM_CONFIG_USE_JENKINS_INSTALL_DEFAULTS}" == "1" ]] ; then
  export ATDM_CONFIG_USE_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE_DEFAULT=1
fi
echo "ATDM_CONFIG_USE_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE_DEFAULT=${ATDM_CONFIG_USE_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE_DEFAULT}"

if  [[ "${ATDM_CONFIG_USE_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR_DEFAULT}" == "" ]] \
  && [[ "${ATDM_CONFIG_USE_JENKINS_INSTALL_DEFAULTS}" == "1" ]] ; then
  export ATDM_CONFIG_USE_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR_DEFAULT=1
fi
echo "ATDM_CONFIG_USE_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR_DEFAULT=${ATDM_CONFIG_USE_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR_DEFAULT}"

if  [[ "${ATDM_CONFIG_USE_MAKE_INSTALL_GROUP_DEFAULT}" == "" ]] \
  && [[ "${ATDM_CONFIG_USE_JENKINS_INSTALL_DEFAULTS}" == "1" ]] ; then
  export ATDM_CONFIG_USE_MAKE_INSTALL_GROUP_DEFAULT=1
fi
echo "ATDM_CONFIG_USE_MAKE_INSTALL_GROUP_DEFAULT=${ATDM_CONFIG_USE_MAKE_INSTALL_GROUP_DEFAULT}"

# Set up the default workspace and base install directory paths

echo

if [[ "${ATDM_CONFIG_WORKSPACE_BASE}" == "" ]] \
  && [[ "${ATDM_CONFIG_WORKSPACE_BASE_DEFAULT}" != "" ]] \
  && [[ "${ATDM_CONFIG_USE_WORKSPACE_BASE_DEFAULT}" == "1" ]] ; then
  export ATDM_CONFIG_WORKSPACE_BASE="${ATDM_CONFIG_WORKSPACE_BASE_DEFAULT}"
fi
echo "ATDM_CONFIG_WORKSPACE_BASE=${ATDM_CONFIG_WORKSPACE_BASE}"

if [[ "${ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE}" == "" ]] \
  && [[ "${ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE_DEFAULT}" != "" ]] \
  && [[ "${ATDM_CONFIG_USE_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE_DEFAULT}" == "1" ]] ; then
  export ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE="${ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE_DEFAULT}"
fi
echo "ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE=${ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE}"

if [ "${ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE}" != "" ] ; then

  echo
  echo "Setting up to use a standard <install-prefix-base>/<date>/<system-build-name>/ install dir!"
  echo

  # Get the <date> (YYYY-MM-DD) as part of directory name.  (NOTE: This <date>
  # should exactly match the CDash 'date=<date>` field in the various PHP
  # pages where this build will end up.
  CDASH_TESTING_DATE=`${WORKSPACE}/Trilinos/cmake/ctest/drivers/trilinos_cdash_build_testing_day.sh`

  if [[ "${ATDM_CONFIG_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR}" == "" ]] \
    && [[ "${ATDM_CONFIG_USE_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR_DEFAULT}" == "1" ]]  ; then
    export ATDM_CONFIG_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR="${ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE}/${CDASH_TESTING_DATE}"
  fi
  echo "ATDM_CONFIG_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR=${ATDM_CONFIG_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR}"

  # Get a unique name for the build that includes the system name, but not the
  # 'Trilinos-atdm-' prefix.  For example, for the build name
  # 'Trilinos-atdm-cee-rhel7_clang-5.0.1_openmpi-1.10.2_serial_static_opt'
  # this pulls out 'cee-rhel7_clang-5.0.1_openmpi-1.10.2_serial_static_opt'.
  # By including the system name in the directory name, this should make the
  # name unique across all platforms allowing a single base directory to be
  # used for builds from many different machines in a mounted drive (for
  # example).
  SYSTEM_AND_BUILD_NAME=`echo $JOB_NAME | sed 's/.*Trilinos-atdm-\(.*\)$/\1/'`

  # Full install dir path <install-prefix-base>/<date>/<system-build-name>
  export ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX="${ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE}/${CDASH_TESTING_DATE}/${SYSTEM_AND_BUILD_NAME}"
  # NOTE: Above, not using
  # ATDM_CONFIG_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR in case that var
  # is empty.

  # Show the full install dir path
  echo "ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX=${ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX}"

fi

if [[ "${ATDM_CONFIG_MAKE_INSTALL_GROUP}" == "" ]] \
  && [[ "${ATDM_CONFIG_MAKE_INSTALL_GROUP_DEFAULT}" != "" ]] \
  && [[ "${ATDM_CONFIG_USE_MAKE_INSTALL_GROUP_DEFAULT}" == "1" ]] ; then
  export ATDM_CONFIG_MAKE_INSTALL_GROUP="${ATDM_CONFIG_MAKE_INSTALL_GROUP_DEFAULT}"
fi
echo "ATDM_CONFIG_MAKE_INSTALL_GROUP=${ATDM_CONFIG_MAKE_INSTALL_GROUP}"

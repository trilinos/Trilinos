#
# A) Load the env
#

source ${WORKSPACE}/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME
echo
module list
echo
echo "cmake in path:"
which cmake
echo
echo "ninja in path:"
which ninja
echo
echo "ATDM config env vars:"
set | grep ATDM_CONFIG_
echo
echo "PATH=$PATH"

#
# B) Set the defaut git repo for Trilinos
#

if [ "${Trilinos_REPOSITORY_LOCATION}" == "" ] ; then
  export Trilinos_REPOSITORY_LOCATION=https://github.com/trilinos/Trilinos.git
fi

unset http_proxy
# NOTE: Above we have to unset http_proxy to allow the second submit to the
# testing-dev.sandia.gov/cdash/ site which the jenkins job sets.  But we can't
# unset https_proxy which is needed for the git operations with
# https://github.com.

#
# C) Setup install-releated stuff
#

echo

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

  # Get a unique name for the build that includes the system name, but not the
  # 'Trilinos-atdm-' prefix.  For example, for the build name
  # 'Trilinos-atdm-cee-rhel6_clang-5.0.1_openmpi-1.10.2_serial_static_opt'
  # this pulls out 'cee-rhel6_clang-5.0.1_openmpi-1.10.2_serial_static_opt'.
  # By including the system name in the directory name, this should make the
  # name unique across all platforms allowing a single base directory to be
  # used for builds from many different machines in a mounted drive (for
  # example).
  SYSTEM_AND_BUILD_NAME=`echo $JOB_NAME | sed 's/.*Trilinos-atdm-\(.*\)$/\1/'`

  # Full install dir path <install-prefix-base>/<date>/<system-build-name>
  export ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX="${ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE}/${CDASH_TESTING_DATE}/${SYSTEM_AND_BUILD_NAME}"

  # Show the full install dir path
  echo "ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX=${ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX}"

fi

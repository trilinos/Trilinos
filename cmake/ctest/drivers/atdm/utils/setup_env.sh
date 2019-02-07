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

if [ "${Trilinos_REPOSITORY_LOCATION}" == "" ] ; then
  export Trilinos_REPOSITORY_LOCATION=https://github.com/trilinos/Trilinos.git
fi

# Echo the env var set by the driver (if desired)
echo "ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE = ${ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE}"

if [ "${ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE}" != "" ] ; then

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
  echo "ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX = ${ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX}"

fi

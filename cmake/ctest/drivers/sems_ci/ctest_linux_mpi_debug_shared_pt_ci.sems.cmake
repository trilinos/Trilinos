
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.sems.cmake")

#
# Set the options specific to this build case
#

SET(BUILD_DIR_NAME MPI_RELEASE_DEBUG_SHARED_PT_OPENMP_CI)
#SET(CTEST_TEST_TIMEOUT 900)

#override the default number of processors to run on and use Ninja
SET( CTEST_CMAKE_GENERATOR Ninja )
SET( CTEST_BUILD_FLAGS "-j16 -k 999999" )
SET( CTEST_PARALLEL_LEVEL "16" )

SET(Trilinos_CTEST_DO_ALL_AT_ONCE TRUE)

SET(Trilinos_ENABLE_SECONDARY_TESTED_CODE OFF)

SET(Trilinos_ENABLE_CONFIGURE_TIMING ON)

SET(Trilinos_BRANCH develop)

SET(EXTRA_EXCLUDE_PACKAGES)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/GCC-4.8.4-OpenMPI-1.10.1-MpiReleaseDebugSharedPtOpenMP.cmake"
  "-DTrilinos_TEST_CATEGORIES=BASIC"
  "-DTrilinos_ENABLE_CONFIGURE_TIMING=ON"
  )
# NOTE: That above must match *exactly* what is listed is listed in
# project-checkin-test-config.py and produced by the checkin-test-sems.sh
# --default-builds=MPI_RELEASE_DEBUG_SHARED_PT_OPENMP build!

SET(CTEST_TEST_TYPE Continuous)

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

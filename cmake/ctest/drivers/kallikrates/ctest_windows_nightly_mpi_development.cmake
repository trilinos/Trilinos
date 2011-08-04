INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.kallikrates.msvc.cmake")

#
# Set the options specific to this build case
#
#SET(CTEST_DO_UPDATES FALSE)
SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME MPI_OPT_DEV)
#SET(CTEST_TEST_TYPE EXPERIMENTAL)
#SET(CTEST_TEST_TIMEOUT 900)

SET( EXTRA_EXCLUDE_PACKAGES TrilinosFramework Amesos Ifpack ML Moertel Stokhos Piro Optika Mesquite FEApp Intrepid)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE:STRING=Experimental"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

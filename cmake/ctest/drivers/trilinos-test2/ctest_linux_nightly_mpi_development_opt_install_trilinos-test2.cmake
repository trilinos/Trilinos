#run command to create a tarball and install from that so we can test against that installation

#hard coding where the installation dir is for now.
#Ideally this should be pulled in from the same variable that is used to set the installation
#dir for the script.
SET(INSTALLATION_DIR "${CMAKE_CURRENT_BINARY_DIR}/../../installation/installed")

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.trilinos-test2.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME MPI_OPT_DEV_INSTALL)
#SET(CTEST_TEST_TIMEOUT 900)
SET(CTEST_TEST_TYPE Nightly)

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE OFF)

SET(EXTRA_EXCLUDE_PACKAGES Mesquite RBGen)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_INSTALLATION_DIR=${INSTALLATION_DIR}"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-DMPI_BASE_DIR:PATH=/home/trilinos"
  "-DTrilinos_ENABLE_DEVELOPMENT_MODE:BOOL=OFF"
  "-DTPL_ENABLE_Pthread:BOOL=ON"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  "-DBoost_INCLUDE_DIRS:FILEPATH=/home/trilinos/include"
  "-DTrilinos_ENABLE_RBGen=OFF"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

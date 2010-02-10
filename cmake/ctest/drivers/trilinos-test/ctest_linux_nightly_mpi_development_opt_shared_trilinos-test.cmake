
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.trilinos-test.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME MPI_OPT_DEV_SHARED)

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)

#disabling Mesquite because of a build error when shared libs is turned on.
SET(EXTRA_EXCLUDE_PACKAGES Mesquite)
SET( Trilinos_EXCLUDE_PACKAGES ${EXTRA_EXCLUDE_PACKAGES} TriKota)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DDART_TESTING_TIMEOUT:STRING=120"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-DBUILD_SHARED_LIBS:BOOL=ON"
  "-DMPI_BASE_DIR:PATH=/home/trilinos/openmpi-1.4"
  "-DTPL_ENABLE_Pthread:BOOL=ON"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  "-DBoost_INCLUDE_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/boost_1_38_0"
  "-DSWIG_EXECUTABLE:FILEPATH=/home/trilinos/tpl/gcc4.1.2/swig-1.3.40/bin/swig"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

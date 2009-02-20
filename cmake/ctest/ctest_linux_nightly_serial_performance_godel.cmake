
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.godel.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE NONE)
SET(BUILD_DIR_NAME SERIAL_PERF)

SET(Trilinos_PACKAGES Teuchos Sacado)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_CXX_FLAGS:STRING=\"-03\""
  "-DCMAKE_C_FLAGS:STRING=\"-03\""
  "-DCMAKE_Fortran_FLAGS:STRING=\"-O5\""
  "-DDART_TESTING_TIMEOUT:STRING=600"
  "-DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE"
  "-DTrilinos_ENABLE_DEPENCENCY_UNIT_TESTS:BOOL=OFF"
  "-DTrilinos_ENABLE_PERFORMANCE_TESTS:BOOL=ON"
  "-DTrilinos_ENABLE_TESTS:BOOL=OFF"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  )

# NOTE: Above, we want to overridde the enabling of all tests and only
# want to do performance tests.

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

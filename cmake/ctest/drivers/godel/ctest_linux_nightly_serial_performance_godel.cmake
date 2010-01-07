
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.godel.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE NONE)
SET(BUILD_DIR_NAME SERIAL_PERF)

SET(Trilinos_PACKAGES Teuchos Sacado)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_CXX_FLAGS:STRING=\"-03 -DBOOST_SP_DISABLE_THREADS\""
  "-DCMAKE_C_FLAGS:STRING=\"-03\""
  "-DCMAKE_Fortran_FLAGS:STRING=\"-O5\""
  "-DDART_TESTING_TIMEOUT:STRING=120"
  "-DTrilinos_ENABLE_DEPENCENCY_UNIT_TESTS:BOOL=OFF"
  "-DTrilinos_ENABLE_PERFORMANCE_TESTS:BOOL=ON"
  "-DTrilinos_ENABLE_TESTS:BOOL=OFF"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  )

# NOTE: Above, we want to overridde the enabling of all tests and only
# want to do performance tests.

# NOTE: We are turning off all optional packages because we don't want
# build errors in other packges to affect these packages and we don't
# need them for performance tests.  The Sacado performance build
# failed on 2009/04/23 because Zoltan failed to build.

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

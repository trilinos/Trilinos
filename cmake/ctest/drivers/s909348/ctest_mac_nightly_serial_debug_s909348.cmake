
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.s909348.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME SERIAL_DEBUG)

# Exclude Sundance because of strange segfault (see bug 4382)
SET( Trilinos_EXCLUDE_PACKAGES Sundance )

SET( EXTRA_CONFIGURE_OPTIONS
  "-DDART_TESTING_TIMEOUT:STRING=600"
  "-DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

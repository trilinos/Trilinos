
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../../TrilinosVersion.cmake")
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.brain.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME "SERIAL_DEBUG_GCC-4.6.1_CXX11_DEV")
SET(COMPILER_VERSION "GCC-4.6.1-C++1X")
SET(Trilinos_TRACK ${Trilinos_TESTING_TRACK})
SET(Trilinos_BRANCH ${Trilinos_REPOSITORY_BRANCH})

# Set up the compiler environment.
SET(brain_GCC_ROOT "$ENV{HOME}/compilers/gcc-4.6.1")

SET(ENV{LD_LIBRARY_PATH} "${brain_GCC_ROOT}/lib64:$ENV{LD_LIBRARY_PATH}")

SET(Trilinos_PACKAGES Tpetra)

SET( EXTRA_CONFIGURE_OPTIONS
  "${EXTRA_CONFIGURE_OPTIONS}"
  "-DTrilinos_ENABLE_CXX11:BOOL=ON"
  "-DCMAKE_C_COMPILER:FILEPATH=${brain_GCC_ROOT}/bin/gcc"
  "-DCMAKE_CXX_COMPILER:FILEPATH=${brain_GCC_ROOT}/bin/g++"
  "-DCMAKE_Fortran_COMPILER:FILEPATH=${brain_GCC_ROOT}/bin/gfortran"
  "-DCMAKE_CXX_FLAGS:STRING=-std=c++0x"
)

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()



INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.exetazo.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME SERIAL_DEBUG_RELEASE_10.0)
#SET(CTEST_TEST_TYPE EXPERIMENTAL)
SET(Trilinos_TRACK "Nightly Release 10.0")
set(CTEST_TEST_TIMEOUT "720")

SET(Trilinos_BRANCH "-r trilinos-release-10-0-branch")

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_VERBOSE_CONFIGURE:BOOL=ON"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTPL_ENABLE_UMFPACK:BOOL=ON"
  "-DTPL_ENABLE_SuperLU:BOOL=ON"
  "-DSuperLU_INCLUDE_DIRS:PATH=/home/jmwille/lib/SuperLU_3.0/SRC"
  "-DSuperLU_LIBRARY_DIRS:PATH=/home/jmwille/lib/SuperLU_3.0"
  "-DSuperLU_LIBRARY_NAMES:STRING=superlu_3.0"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  "-DBoost_INCLUDE_DIRS:FILEPATH=/home/jmwille/boost_1_35_0"
  "-DTPL_ENABLE_Oski:BOOL=ON"
  "-DTPL_ENABLE_y12m:BOOL=ON"
  "-Dy12m_LIBRARY_DIRS:FILEPATH=/home/jmwille/install/y12m/lib"
)

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

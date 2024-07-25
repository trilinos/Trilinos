  
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for trilinos-test2 using gcc
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME
  SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )

  SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  
  SET_DEFAULT( CTEST_BUILD_FLAGS "-j10 -i" )

  SET_DEFAULT( CTEST_PARALLEL_LEVEL "10" )

  SET_DEFAULT( Trilinos_ENABLE_SECONDARY_TESTED_CODE OFF)

  SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS )
  # No options to set!  When the SEMS env is loaded correctly, the compilers,
  # MPI, and the TPLs will be found automatically!

  SET_DEFAULT(COMPILER_VERSION "GCC-4.8.4")

  TRILINOS_CTEST_DRIVER()

ENDMACRO()

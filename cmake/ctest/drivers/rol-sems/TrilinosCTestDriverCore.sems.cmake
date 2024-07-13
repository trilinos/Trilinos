  
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for trilinos-test2 using gcc
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME
  SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )

  SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  
  SET_DEFAULT( CTEST_BUILD_FLAGS "-j20 -i" )

  SET_DEFAULT( CTEST_PARALLEL_LEVEL "20" )

  SET_DEFAULT( Trilinos_ENABLE_SECONDARY_TESTED_CODE OFF)

  SET_DEFAULT( Trilinos_PACKAGES "ROL" )

  SET_DEFAULT( CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES "OFF" )

  SET_DEFAULT( Trilinos_TRACK Specialized )

  #SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS )
  # No options to set!  When the SEMS env is loaded correctly, the compilers,
  # MPI, and the TPLs will be found automatically!

  SET_DEFAULT( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/MpiReleaseDebugSharedPtSettings.cmake,cmake/std/BasicCiTestingSettings.cmake,cmake/std/sems/SEMSDevEnv.cmake"
  "-DTrilinos_TEST_CATEGORIES=BASIC"
  "-DTrilinos_ENABLE_CONFIGURE_TIMING=ON"
  )


  SET_DEFAULT(COMPILER_VERSION "GCC-4.8.4")

  TRILINOS_CTEST_DRIVER()

ENDMACRO()

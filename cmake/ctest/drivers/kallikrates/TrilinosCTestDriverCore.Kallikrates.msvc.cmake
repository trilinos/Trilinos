  
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for godel using gcc
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME
  SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )

  SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  # convert CVS_EXE path to a cmake style path and store
  # in CVS_EXECUTABLE. CVS_EXE is an environment variable
  # that should be set before running this script.
  FILE(TO_CMAKE_PATH "$ENV{CVS_EXE}" CVS_EXECUTABLE)
  SET( CTEST_CMAKE_GENERATOR "NMake Makefiles")
  SET( CTEST_BUILD_FLAGS " -i" )
  SET_DEFAULT( Trilinos_ENABLE_SECONDARY_STABLE_CODE ON )
  SET(COMPILER_VERSION "MSVC9")

  SET_DEFAULT( Trilinos_EXCLUDE_PACKAGES ${EXTRA_EXCLUDE_PACKAGES} )

  SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
    "-DCLAPACK_DIR=C:/trilinos_projects/clapack_build"
    "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
    "-DBoost_INCLUDE_DIRS=C:/trilinos_projects/boost_1_40_0"
    "-DTPL_ENABLE_Netcdf=OFF"
    )

  IF (COMM_TYPE STREQUAL MPI)

    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DTPL_ENABLE_MPI:BOOL=ON"
      )

  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()

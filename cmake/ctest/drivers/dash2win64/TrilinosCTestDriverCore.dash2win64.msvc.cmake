
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for Kitware Windows machine dash2win64
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
  SET_DEFAULT(Trilinos_EXCLUDE_PACKAGES PyTrilinos TriKota Optika)
  SET(COMPILER_VERSION "MSVC9")
  SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
    "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
    "-DDART_TESTING_TIMEOUT:STRING=120"
    "-DBoost_INCLUDE_DIRS:PATH=C:/Dashboards/Support/boost_1_39_0"
    "-DTrilinos_ENABLE_TriKota:BOOL=OFF"
    )

  TRILINOS_CTEST_DRIVER()

ENDMACRO()

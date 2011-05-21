  
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for godel using gcc
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME
  SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )

  SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  
  SET( CTEST_BUILD_FLAGS "-j8 -i" )
  SET( CTEST_PARALLEL_LEVEL 8 )

  #SET( CTEST_MEMORYCHECK_COMMAND /usr/bin/valgrind )
  #SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS )

  # We don't have TPLs for PyTrilinos, Optika, and TriKota
  # There is some strange problem with TrilinosFramework causing an strange error.
  SET_DEFAULT( Trilinos_EXCLUDE_PACKAGES TrilinosFramework PyTrilinos TriKota Optika)

  # Enable LIME since it is EX currently
  SET(Trilinos_ADDITIONAL_PACKAGES LIME CASLBOA)
  
  SET_DEFAULT(COMPILER_VERSION "ICPC-12.191")
  
  IF (COMM_TYPE STREQUAL MPI)

    MESSAGE(FATAL_ERROR "Error, Intel build does not support MPI yet!")
  
  ELSE()
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      "-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH=${CTEST_SCRIPT_DIRECTORY}/intel-12.191-options.cmake"
      "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
      "-DTrilinos_ENABLE_TriKota:BOOL=OFF"
    )

  #   "-DCMAKE_CXX_FLAGS:STRING=-diag-disable 597"
  #   "-DCMAKE_LIBRARY_PATH:PATH=/usr/lib64"
  
  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()

  
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
  # 
  SET_DEFAULT( Trilinos_EXCLUDE_PACKAGES
     PyTrilinos TriKota Optika  # We don't have TPLs for these
     Sundance Stokhos STK  # Currently have failures and nor currently needed by CASL
     )
  
  SET_DEFAULT(COMPILER_VERSION "GCC-4.5.1")

  SET_DEFAULT(Trilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE Nightly)
  
  IF (COMM_TYPE STREQUAL MPI)
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      "-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH=${CTEST_SCRIPT_DIRECTORY}/gcc-4.5.1-mpi-ss-options.cmake"
      "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
      "-DTPL_ENABLE_MPI:BOOL=ON"
      "-DTrilinos_ENABLE_TriKota:BOOL=OFF"
    )
  
  ELSE()
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      "-DTrilinos_CONFIGURE_OPTIONS_FILE:FILEPATH=${CTEST_SCRIPT_DIRECTORY}/gcc-4.5.1-serial-ss-options.cmake"
      "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
      "-DTrilinos_ENABLE_TriKota:BOOL=OFF"
    )

  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()

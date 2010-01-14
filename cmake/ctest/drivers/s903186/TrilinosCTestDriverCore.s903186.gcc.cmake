  
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for s903186 using gcc
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME

  SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )

  SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  
  SET( CTEST_BUILD_FLAGS "-j4 -i" )

  SET_DEFAULT( Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)
  
  # Only turn on PyTrilinos for shared libraries
  SET_DEFAULT(Trilinos_EXCLUDE_PACKAGES ${EXTRA_EXCLUDE_PACKAGES} PyTrilinos TriKota Optika)
  
  SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
    "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
    "-DTrilinos_ENABLE_DEPENCENCY_UNIT_TESTS:BOOL=OFF"
    "-DBoost_INCLUDE_DIRS:FILEPATH=/Users/jmwille/install/boost_1_41_0"
    "-DTrilinos_ENABLE_TriKota:BOOL=OFF"
    "-DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE"
    )

  SET_DEFAULT(COMPILER_VERSION "GCC-4.4.2")
  
  IF (COMM_TYPE STREQUAL MPI)
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DTPL_ENABLE_MPI:BOOL=ON"
      "-DMPI_BASE_DIR:PATH=/Users/jmwille/install"
      )
  
  ELSE()
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DCMAKE_CXX_COMPILER:FILEPATH=/Users/jmwille/install/bin/g++"
      "-DCMAKE_C_COMPILER:FILEPATH=/Users/jmwille/install/bin/gcc"
      "-DCMAKE_Fortran_COMPILER:FILEPATH=/Users/jmwille/install/bin/gfortran"
      )
  
  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()

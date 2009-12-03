
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for godel using gcc
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME
  SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )

  SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  
  SET( CTEST_BUILD_FLAGS "-j4 -i" )

  SET( CTEST_COVERAGE_COMMAND /usr/bin/gcov )
  SET( CTEST_MEMORYCHECK_COMMAND /usr/bin/valgrind )

  SET_DEFAULT( Trilinos_ENABLE_SECONDARY_STABLE_CODE OFF )

  # The CheckinTest unit test fails for some reason on gabriel in nighlty
  # mode.  However, it passes when I run it it locally and it gets run nightly
  # on other machines so I will disable this for now.  If someone else reports
  # a problem then I will look into this further.
  SET_DEFAULT( Trilinos_EXCLUDE_PACKAGES TrilinosFramework )
  
  SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
    "-DCMAKE_BUILD_TYPE:STRING=RELEASE"
    "-DDART_TESTING_TIMEOUT:STRING=120"
    "-DCOVERAGE_COMMAND:FILEPATH=/usr/bin/gcov"
    "-DMEMORYCHECK_COMMAND:FILEPATH=/usr/bin/valgrind"
    "-DZoltan_ENABLE_TESTS:BOOL=OFF"
    )

  IF (BUILD_TYPE STREQUAL "DEBUG")
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DTrilinos_ENABLE_DEBUG:BOOL=ON"
      )
  ENDIF()
  
  IF (COMM_TYPE STREQUAL MPI)
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DTPL_ENABLE_MPI:BOOL=ON"
      )

    SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS
        "--trace-children=yes --gen-suppressions=all --suppressions=${CTEST_SCRIPT_DIRECTORY}/valgrind_suppressions_godel_openmpi_1.2.7.txt  ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}" )

    SET_DEFAULT(COMPILER_VERSION "GCC-3.4.6")
  
  ELSE()
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      #"-DCMAKE_CXX_COMPILER:FILEPATH=/usr/local/gcc-4.2.0/bin/g++"
      #"-DCMAKE_C_COMPILER:FILEPATH=/usr/local/gcc-4.2.0/bin/gcc"
      #"-DCMAKE_Fortran_COMPILER:FILEPATH=/usr/local/gcc-4.2.0/bin/gfortran"
      )

    SET_DEFAULT(COMPILER_VERSION "GCC-3.4.6")
    #SET_DEFAULT(COMPILER_VERSION "GCC-4.2.0")
  
  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()

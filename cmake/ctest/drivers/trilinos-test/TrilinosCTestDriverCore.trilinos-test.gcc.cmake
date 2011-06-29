  
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for trilinos-test using gcc
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME
  SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )

  SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  
  SET_DEFAULT( CTEST_BUILD_FLAGS "-j8 -i" )

  SET_DEFAULT( CTEST_PARALLEL_LEVEL "8" )

  SET( CTEST_COVERAGE_COMMAND /usr/bin/gcov )
  SET( CTEST_MEMORYCHECK_COMMAND /usr/bin/valgrind )

  SET_DEFAULT( Trilinos_ENABLE_SECONDARY_STABLE_CODE OFF)

  # Only turn on PyTrilinos for shared libraries
  SET_DEFAULT( Trilinos_EXCLUDE_PACKAGES ${EXTRA_EXCLUDE_PACKAGES} PyTrilinos TriKota)
  
  SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
    "-DCOVERAGE_COMMAND:FILEPATH=/usr/bin/gcov"
    "-DMEMORYCHECK_COMMAND:FILEPATH=/usr/bin/valgrind"
    "-DTPL_ENABLE_Netcdf:BOOL=ON"
    "-DNetcdf_LIBRARY_DIRS=/home/trilinos/tpl/gcc4.1.2/netcdf_4.0/lib"
    "-DNetcdf_INCLUDE_DIRS=/home/trilinos/tpl/gcc4.1.2/netcdf_4.0/include"
    "-DTPL_ENABLE_SuperLU:BOOL=ON"
    "-DSuperLU_INCLUDE_DIRS:PATH=/home/trilinos/tpl/gcc4.1.2/SuperLU_4.1/SRC"
    "-DSuperLU_LIBRARY_DIRS:PATH=/home/trilinos/tpl/gcc4.1.2/SuperLU_4.1/lib"
    "-DSuperLU_LIBRARY_NAMES:STRING=superlu_4.1"
    "-DTrilinos_ENABLE_TriKota:BOOL=OFF"
    "-DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE"
    "-DMesquite_ENABLE_TESTS:BOOL=ON"
    "-DCPPUNIT_LIBRARY:STRING=/home/trilinos/tpl/gcc4.1.2/cppunit-1.12.1/lib/libcppunit.a"
    "-DCPPUNIT_INCLUDES:STRING=/home/trilinos/tpl/gcc4.1.2/cppunit-1.12.1/include"
    )

  IF (BUILD_TYPE STREQUAL "DEBUG")
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DTrilinos_ENABLE_DEBUG:BOOL=ON"
      )
  ENDIF()

  SET_DEFAULT(COMPILER_VERSION "GCC-4.1.2")
  
  IF (COMM_TYPE STREQUAL MPI)
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DTPL_ENABLE_MPI:BOOL=ON"
      )
#      "-DMPI_BASE_DIR:PATH=/home/trilinos/openmpi-1.4"
#      )

    SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS
        "--trace-children=yes --gen-suppressions=all --suppressions=${CTEST_SCRIPT_DIRECTORY}/valgrind_suppressions_trilinos-test_openmpi_1.2.7.txt ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}" )
  
  ELSE()
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++"
      "-DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc"
      "-DCMAKE_Fortran_COMPILER:FILEPATH=/usr/bin/gfortran"
      )

    SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS
        "--trace-children=yes --gen-suppressions=all --suppressions=${CTEST_SCRIPT_DIRECTORY}/valgrind_suppressions_trilinos-test_gcc-4.1.2.txt ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}" )
  
  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()

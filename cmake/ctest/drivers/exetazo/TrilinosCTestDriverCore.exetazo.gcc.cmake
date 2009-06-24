  
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for exetazo using gcc
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME
  SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )

  SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  
  SET( CTEST_CVS_COMMAND_ARGS "cvs -q -z3" )
  
  SET( CTEST_BUILD_FLAGS "-j6 -i" )

  SET( CTEST_COVERAGE_COMMAND /usr/bin/gcov )
  SET( CTEST_MEMORYCHECK_COMMAND /usr/bin/valgrind )

  SET_DEFAULT( Trilinos_ENABLE_SECONDARY_STABLE_CODE OFF)

  # Only turn on PyTrilinos for shared libraries
  SET_DEFAULT( Trilinos_EXCLUDE_PACKAGES PyTrilinos )
  
  SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
    "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
    "-DCOVERAGE_COMMAND:FILEPATH=/usr/bin/gcov"
    "-DMEMORYCHECK_COMMAND:FILEPATH=/usr/bin/valgrind"
    )

  SET_DEFAULT(COMPILER_VERSION "GCC-4.1.2")
  
  IF (COMM_TYPE STREQUAL MPI)
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DTPL_ENABLE_MPI:BOOL=ON"
      "-D MPIEXEC_MAX_NUMPROCS:STRING=2"
      )
#      "-DMPI_BASE_DIR:PATH=/usr/lib64/openmpi/1.2.7-gcc"
#      )

    SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS
        "--trace-children=yes --gen-suppressions=all --suppressions=${CTEST_SCRIPT_DIRECTORY}/valgrind_suppressions_exetazo_openmpi_1.2.7.txt ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}" )
  
  ELSE()
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++"
      "-DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc"
      "-DCMAKE_Fortran_COMPILER:FILEPATH=/usr/bin/gfortran"
      )

    SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS
        "--trace-children=yes --gen-suppressions=all --suppressions=${CTEST_SCRIPT_DIRECTORY}/valgrind_suppressions_exetazo_gcc-4.1.2.txt ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}" )
  
  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()

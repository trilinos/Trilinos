  
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for brain.sandia.gov using gcc 4.1.2
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.brain.common.cmake")

  SET_DEFAULT(COMPILER_VERSION "GCC-4.1.2")

  SET( CTEST_MEMORYCHECK_COMMAND /usr/bin/valgrind )

  SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
    ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
    "-DMEMORYCHECK_COMMAND:FILEPATH=/usr/bin/valgrind"
    )
      
  IF (COMM_TYPE STREQUAL MPI)
  
    # Set MPI_BASE_DIR to use the wrapper compilers.
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DTPL_ENABLE_MPI:BOOL=ON"
      "-DMPI_BASE_DIR=/usr/local"
      )

    SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS
        "--trace-children=yes --gen-suppressions=all --suppressions=${CTEST_SCRIPT_DIRECTORY}/valgrind_suppressions_brain_openmpi_1.2.7.txt ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}" )

  ELSE()
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      )

    SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS
        "--trace-children=yes --gen-suppressions=all --suppressions=${CTEST_SCRIPT_DIRECTORY}/valgrind_suppressions_brain_gcc-4.1.2.txt ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}" )
  
  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()

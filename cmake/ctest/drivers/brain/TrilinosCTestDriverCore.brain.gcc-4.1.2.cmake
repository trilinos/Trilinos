  
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for brain.sandia.gov using gcc 4.1.2
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.brain.common.cmake")

  SET_DEFAULT(COMPILER_VERSION "GCC-4.1.2")
  
  IF (COMM_TYPE STREQUAL MPI)
  
    # Set MPI_BASE_DIR to use the wrapper compilers.
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DTPL_ENABLE_MPI:BOOL=ON"
      "-DMPI_BASE_DIR=/usr/local"
      )
  
  ELSE()
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      )
  
  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()

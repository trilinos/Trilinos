  
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for brain.sandia.gov using gcc 4.6.3
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.brain.common.cmake")

  # Set up the compiler environment.
  SET(brain_GCC_ROOT "$ENV{HOME}/compilers/gcc/4.6.3")
  SET(ENV{LD_LIBRARY_PATH} "${brain_GCC_ROOT}/lib64:$ENV{LD_LIBRARY_PATH}")
  
  SET_DEFAULT(COMPILER_VERSION "GCC-4.6.3")
  
  IF (COMM_TYPE STREQUAL MPI)
  
    # Set MPI_BASE_DIR to use the wrapper compilers.
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DTPL_ENABLE_MPI:BOOL=ON"
      "-DMPI_BASE_DIR=$ENV{HOME}/tpls/gcc/4.6.3/openmpi/1.4.3"
      )
  
  ELSE()
  
    # Explicitly set the compiler variables.
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DCMAKE_C_COMPILER:FILEPATH=${brain_GCC_ROOT}/bin/gcc"
      "-DCMAKE_CXX_COMPILER:FILEPATH=${brain_GCC_ROOT}/bin/g++"
      "-DCMAKE_Fortran_COMPILER:FILEPATH=${brain_GCC_ROOT}/bin/gfortran"
      )
  
  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()

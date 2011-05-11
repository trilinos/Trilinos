  
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

  # Turn off SS packages we don't have TPLs for
  SET_DEFAULT( Trilinos_EXCLUDE_PACKAGES PyTrilinos TriKota Optika)
  
  SET_DEFAULT(COMPILER_VERSION "ICPC-12.0.2")
  
  IF (COMM_TYPE STREQUAL MPI)

    MESSAGE(FATAL_ERROR "Error, Intel build does not support MPI yet!")
  
  ELSE()
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
      "-DTrilinos_ENABLE_TriKota:BOOL=OFF"
      "-DCMAKE_C_COMPILER:FILEPATH=/opt/intel/Compiler/composerxe-2011.2.137/bin/intel64/icc"
      "-DCMAKE_CXX_COMPILER:FILEPATH=/opt/intel/Compiler/composerxe-2011.2.137/bin/intel64/icpc"
      "-DCMAKE_Fortran_COMPILER:FILEPATH=/opt/intel/Compiler/composerxe-2011.2.137/bin/intel64/ifort"
      "-DTPL_BLAS_LIBRARIES:STRING=-L/opt/intel/Compiler/composerxe-2011.2.137/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_blas95_lp64 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_core -lmkl_sequential"
      "-DTPL_LAPACK_LIBRARIES:STRING=-L/opt/intel/Compiler/composerxe-2011.2.137/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_blas95_lp64 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_core -lmkl_sequential"
      "-DTPL_Boost_INCLUDE_DIRS:FILEPATH=/opt/tpls_src/boost_1_46_1"
    )

  #   "-DCMAKE_CXX_FLAGS:STRING=-diag-disable 597"
  #   "-DCMAKE_LIBRARY_PATH:PATH=/usr/lib64"
  
  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()



SET(ARGS "--filedir=${CMAKE_CURRENT_BINARY_DIR}/")

IF (SITE STREQUAL "gabriel.sandia.gov")
  SET(ARGS
    "${ARGS} --not-unit-test=VbrMatrix_int_ComplexFloat_FullMatrixComplex_UnitTest")
ENDIF()

IF(Kokkos_ENABLE_Cuda)
  #
  # Add subpackage include path dependencies for CUDA_COMPILE().
  CUDA_INCLUDE_DIRECTORIES( "." )
  CUDA_INCLUDE_DIRECTORIES( ${Tpetra_INCLUDE_DIRS} )
  CUDA_COMPILE(VM_UT_OBJECTS_CUDA VbrMatrix_UnitTests_CU.cu)
ENDIF()

TRIBITS_ADD_EXECUTABLE(
  VbrMatrix_UnitTests
  SOURCES
    VbrMatrix_UnitTests
    ${VM_UT_OBJECTS_CUDA}
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  )

TRIBITS_ADD_TEST(
  VbrMatrix_UnitTests
  NAME VbrMatrix_UnitTests1
  ARGS ${ARGS}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  VbrMatrix_UnitTests
  NAME VbrMatrix_UnitTests2
  ARGS ${ARGS}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  VbrMatrix_UnitTests
  NAME VbrMatrix_UnitTests3
  ARGS ${ARGS}
  COMM mpi
  NUM_MPI_PROCS 3
  STANDARD_PASS_OUTPUT
  )


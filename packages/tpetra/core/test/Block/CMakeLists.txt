TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExpBlockMultiVector
  SOURCES
    ExpBlockMultiVector
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExpBlockMultiVector2
  SOURCES
    ExpBlockMultiVector2
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExpBlockCrsMatrix
  SOURCES
    ExpBlockCrsMatrix
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(TpetraCore_ExpBlockCrsMatrix_cp
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE_FILES blockA.mm blockA-complex.mm
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExpBlockMap
  SOURCES
    ExpBlockMap
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

#
# As a fix for #506, this test has its own "main" function, to control
# when Kokkos::initialize and Kokkos::finalize get called.  Thus, it
# does not and should not link against ${TEUCHOS_STD_UNIT_TEST_MAIN}.
#
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExpBlockView
  SOURCES
    ExpBlockView
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExpBlockOps
  SOURCES
    ExpBlockOps
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

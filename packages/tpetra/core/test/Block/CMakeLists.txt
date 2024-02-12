TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockMultiVector
  SOURCES
    BlockMultiVector.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockMultiVector2
  SOURCES
    BlockMultiVector2.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockCrsMatrix
  SOURCES
    BlockCrsMatrix.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Issue11931
  SOURCES
  Issue11931.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )



TRIBITS_COPY_FILES_TO_BINARY_DIR(TpetraCore_BlockCrsMatrix_cp
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE_FILES blockA.mm blockA-complex.mm
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockMap
  SOURCES
    BlockMap.cpp
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
  BlockView
  SOURCES
    BlockView.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockOps
  SOURCES
    BlockOps.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )




SET(ARGS "--filedir=${CMAKE_CURRENT_BINARY_DIR}/")

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_UnitTests0
  SOURCES
    CrsGraph_UnitTests0
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  ${MAXNP}
  COMM serial mpi
  FAIL_REGULAR_EXPRESSION "  FAILED  "
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_UnitTests1
  SOURCES
    CrsGraph_UnitTests1
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  ${MAXNP}
  COMM serial mpi
  FAIL_REGULAR_EXPRESSION "  FAILED  "
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_ReindexColumns
  SOURCES
    CrsGraph_ReindexColumns
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS
  COMM serial mpi
  )

# TRIBITS_COPY_FILES_TO_BINARY_DIR(TpetraCrsGraphCopyFiles
#   SOURCE_FILES west0067.rua mhd1280b.cua
#   EXEDEPS CrsGraph_UnitTests
#   )

# This test needs at least 3 MPI processes to be meaningful.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Issue601
  SOURCES
    CrsGraph_Issue601.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM mpi
  NUM_MPI_PROCS 3-8
  STANDARD_PASS_OUTPUT
  )

# This test needs exactly 2 MPI processes to work.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_insertGlobalIndicesFiltered
  SOURCES
    CrsGraph_insertGlobalIndicesFiltered.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

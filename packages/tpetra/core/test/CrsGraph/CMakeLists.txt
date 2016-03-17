

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

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_findRelOffset
  SOURCES
    CrsGraph_findRelOffset
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

# TRIBITS_COPY_FILES_TO_BINARY_DIR(TpetraCrsGraphCopyFiles
#   SOURCE_FILES west0067.rua mhd1280b.cua
#   EXEDEPS CrsGraph_UnitTests
#   )

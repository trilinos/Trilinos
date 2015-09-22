

SET(ARGS "--filedir=${CMAKE_CURRENT_BINARY_DIR}/")

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_UnitTests
  SOURCES
    CrsGraph_UnitTests
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

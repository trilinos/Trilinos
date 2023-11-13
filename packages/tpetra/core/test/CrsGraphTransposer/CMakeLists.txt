
# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   RowMatrixTransposer_test
#   SOURCES
#   main.cpp
#   )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraphTransposer_UnitTests
  SOURCES
    CrsGraphTransposer_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   CrsMatrix_transpose_sortedRows
#   SOURCES sorted.cpp
#   COMM serial mpi
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_COPY_FILES_TO_BINARY_DIR(TransposeGraphFileCopies
  SOURCE_FILES L.mtx LT.mtx L+LT.mtx
  EXEDEPS CrsGraphTransposer_UnitTests
  )

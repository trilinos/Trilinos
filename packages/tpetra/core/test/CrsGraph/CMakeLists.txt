

SET(ARGS "--filedir=${CMAKE_CURRENT_BINARY_DIR}/")

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_UnitTests0
  SOURCES
    CrsGraph_UnitTests0.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  ${MAXNP}
  COMM serial mpi
  FAIL_REGULAR_EXPRESSION "  FAILED  "
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_UnitTests1
  SOURCES
    CrsGraph_UnitTests1.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  ${MAXNP}
  COMM serial mpi
  FAIL_REGULAR_EXPRESSION "  FAILED  "
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_UnitTests_Swap
  SOURCES
  CrsGraph_UnitTests_Swap.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    ARGS "${ARGS} --show-test-details=ALL --output-show-proc-rank"
  ${MAXNP}
  COMM mpi
  NUM_MPI_PROCS 2
  FAIL_REGULAR_EXPRESSION "  FAILED  "
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_ReindexColumns
  SOURCES
    CrsGraph_ReindexColumns.cpp
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

# I want this test to _build_ in either MPI or non-MPI ("serial"),
# but I only want it to _run_ in an MPI build with exactly 1 MPI
# process.
TRIBITS_ADD_EXECUTABLE(
  CrsGraph_PackUnpack
  SOURCES
    CrsGraph_PackUnpack.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
)

TRIBITS_ADD_TEST(
  CrsGraph_PackUnpack
  NAME CrsGraph_PackUnpack_MPI_1
  ARGS ""
  COMM mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_getNumDiags
  SOURCES
    CrsGraph_getNumDiags.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DetermineLocalTriangularStructure
  SOURCES
    DetermineLocalTriangularStructure.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_UnpackIntoStaticGraph
  SOURCES
    CrsGraph_UnpackIntoStaticGraph.cpp
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_StaticImportExport
  SOURCES
    CrsGraph_StaticImportExport.cpp
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_UnpackMerge
  SOURCES
    UnpackMerge.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_DetectNonlocals
  SOURCES
  DetectNonlocals.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    ARGS "--show-test-details=ALL --output-show-proc-rank"
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

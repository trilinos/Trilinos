

SET(ARGS "--filedir=${CMAKE_CURRENT_BINARY_DIR}/")

IF (SITE STREQUAL "gabriel.sandia.gov")
  SET(ARGS
    "${ARGS} --not-unit-test=CrsMatrix_int_ComplexFloat_FullMatrixComplex_UnitTest")
ENDIF()

JOIN(allnodes   "|" FALSE ${Tpetra_ETI_NODES}  )
TRIBITS_ETI_TYPE_EXPANSION(conversions "SOUT=int" "SIN=double" "LO=int" "GO=int" "N=${allnodes}")
TRIBITS_ADD_ETI_INSTANTIATIONS(Tpetra ${conversions})

#
# We break up CrsMatrix's unit tests into several files, to reduce
# compilation time for each file and improve build parallelism.
#

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnitTests
  SOURCES
    CrsMatrix_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnitTests2
  SOURCES
    CrsMatrix_UnitTests2
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnitTests3
  SOURCES
    CrsMatrix_UnitTests3
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnitTests4
  SOURCES
    CrsMatrix_UnitTests4
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_NonlocalAfterResume
  SOURCES
    CrsMatrix_NonlocalAfterResume
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_NodeConversion
  SOURCES
    CrsMatrix_NodeConversion
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_TriSolve
  SOURCES
    CrsMatrix_TriSolve
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_LeftRightScale
  SOURCES
    CrsMatrix_LeftRightScale
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_WithGraph
  SOURCES
    CrsMatrix_WithGraph
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_ReplaceDomainMapAndImporter
  SOURCES
    CrsMatrix_ReplaceDomainMapAndImporter
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

IF (${PACKAGE_NAME}_ENABLE_CLASSIC_VBR)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    CrsMatrix_BlockExtraction
    SOURCES
    CrsMatrix_BlockExtraction
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    # ARGS ${ARGS}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF ()

TRIBITS_COPY_FILES_TO_BINARY_DIR(CrsMatrixCopyFiles1
  SOURCE_FILES west0067.rua mhd1280b.cua
  EXEDEPS CrsMatrix_UnitTests
  )

ASSERT_DEFINED(Anasazi_SOURCE_DIR)
TRIBITS_COPY_FILES_TO_BINARY_DIR(CrsTests_CopyFiles2
  SOURCE_DIR ${Anasazi_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk14.hb bcsstk15.hb
  )

# This test only makes sense for > 1 processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_NonlocalSumInto
  SOURCES
    CrsMatrix_NonlocalSumInto
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS
  COMM mpi
  NUM_MPI_PROCS 2-4
  STANDARD_PASS_OUTPUT
  )

# This test only makes sense for > 1 processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_NonlocalSumInto_Ignore
  SOURCES
    CrsMatrix_NonlocalSumInto_Ignore
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS
  COMM mpi
  NUM_MPI_PROCS 2-4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_gaussSeidel
  SOURCES
    CrsMatrix_gaussSeidel
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )

# This test only makes sense for exactly 2 MPI processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_Bug5978
  SOURCES
    bug5978
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

# First test for Bug 6069.  This test only makes sense for exactly 3
# MPI processes.  Run under Valgrind for best results.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_Bug6069_1
  SOURCES
    Bug6069_1
  ARGS
  COMM mpi
  NUM_MPI_PROCS 3
  STANDARD_PASS_OUTPUT
  )

# Second test for Bug 6069.  This test only makes sense for exactly 2
# MPI processes.  Run under Valgrind for best results.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_Bug6069_2
  SOURCES
    Bug6069_2
  ARGS
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

# Test for Bug 6171.  This test only makes sense for exactly 2 MPI
# processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_Bug6171
  SOURCES
    Bug6171
  ARGS ""
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_ReplaceLocalValues
  SOURCES
    CrsMatrix_ReplaceLocalValues
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_MultipleFillCompletes
  SOURCES
    CrsMatrix_MultipleFillCompletes
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_ReindexColumns
  SOURCES
    CrsMatrix_ReindexColumns
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )


SET(TIMING_INSTALLS "")

INSTALL(TARGETS ${TIMING_INSTALLS}
        RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")

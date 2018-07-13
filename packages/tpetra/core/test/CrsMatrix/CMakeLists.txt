ADD_SUBDIRECTORIES(Regression)

SET(ARGS "--filedir=${CMAKE_CURRENT_BINARY_DIR}/")

#
# We break up CrsMatrix's unit tests into several files, to reduce
# compilation time for each file and improve build parallelism.
#

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnitTests
  SOURCES
    CrsMatrix_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
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
  CrsMatrix_LeftRightScale
  SOURCES
    CrsMatrix_LeftRightScale
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# We split the CrsMatrix_WithGraph test by execution space.
# This speeds up the build.

ASSERT_DEFINED (Tpetra_INST_CUDA)
IF (Tpetra_INST_CUDA)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    CrsMatrix_WithGraph_Cuda
    SOURCES
      CrsMatrix_WithGraph_Cuda
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    ARGS ${ARGS}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF ()

ASSERT_DEFINED (Tpetra_INST_OPENMP)
IF (Tpetra_INST_OPENMP)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    CrsMatrix_WithGraph_OpenMP
    SOURCES
      CrsMatrix_WithGraph_OpenMP
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    ARGS ${ARGS}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF ()

ASSERT_DEFINED (Tpetra_INST_SERIAL)
IF (Tpetra_INST_SERIAL)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    CrsMatrix_WithGraph_Serial
    SOURCES
      CrsMatrix_WithGraph_Serial
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    ARGS ${ARGS}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF ()

ASSERT_DEFINED (Tpetra_INST_PTHREAD)
IF (Tpetra_INST_PTHREAD)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    CrsMatrix_WithGraph_Threads
    SOURCES
      CrsMatrix_WithGraph_Threads
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    ARGS ${ARGS}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF ()


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_ReplaceDomainMapAndImporter
  SOURCES
    CrsMatrix_ReplaceDomainMapAndImporter
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(CrsMatrixCopyFiles1
  SOURCE_FILES west0067.rua mhd1280b.cua
  EXEDEPS CrsMatrix_UnitTests
  )

ASSERT_DEFINED(Anasazi_SOURCE_DIR)
TRIBITS_COPY_FILES_TO_BINARY_DIR(CrsTests_CopyFiles2
  SOURCE_DIR ${Anasazi_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk14.hb
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
  Bug5978
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
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
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
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
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
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS
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

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_TransformValues
  SOURCES
    CrsMatrix_TransformValues
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_GetRowCopy
  SOURCES
    CrsMatrix_GetRowCopy
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# I want this test to _build_ in either MPI or non-MPI ("serial"),
# but I only want it to _run_ in an MPI build with exactly 1 MPI
# process.
TRIBITS_ADD_EXECUTABLE(
  CrsMatrix_PackUnpack
  SOURCES
    CrsMatrix_PackUnpack.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
)

TRIBITS_ADD_TEST(
  CrsMatrix_PackUnpack
  NAME CrsMatrix_PackUnpack_MPI_1
  ARGS ""
  COMM mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  assembleElement
  SOURCES
    assembleElement
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Equilibration
  SOURCES 
    Equilibration.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

SET(TIMING_INSTALLS "")

INSTALL(TARGETS ${TIMING_INSTALLS}
        RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")

ADD_SUBDIRECTORIES(Regression)

SET(ARGS "--filedir=${CMAKE_CURRENT_BINARY_DIR}/")

#
# We break up CrsMatrix's unit tests into several files, to reduce
# compilation time for each file and improve build parallelism.
#

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Bug10008
  SOURCES
    Bug10008.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnitTests
  SOURCES
    CrsMatrix_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnitTests2
  SOURCES
    CrsMatrix_UnitTests2.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnitTests3
  SOURCES
    CrsMatrix_UnitTests3.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnitTests4
  SOURCES
    CrsMatrix_UnitTests4.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnitTests_LocalAccessors
  SOURCES
    CrsMatrix_UnitTests_LocalAccessors.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnitTests_Swap
  SOURCES
    CrsMatrix_UnitTests_Swap.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_NonlocalAfterResume
  SOURCES
    CrsMatrix_NonlocalAfterResume.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_NodeConversion
  SOURCES
    CrsMatrix_NodeConversion.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_LeftRightScale
  SOURCES
    CrsMatrix_LeftRightScale.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_2DRandomDist
  SOURCES
    CrsMatrix_2DRandomDist.cpp
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

# We split the CrsMatrix_WithGraph test by execution space.
# This speeds up the build.

ASSERT_DEFINED (Tpetra_INST_SYCL)
IF (Tpetra_INST_SYCL)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
	  CrsMatrix_WithGraph_SYCL
    SOURCES
      CrsMatrix_WithGraph_SYCL.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    ARGS ${ARGS}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF ()

ASSERT_DEFINED (Tpetra_INST_CUDA)
IF (Tpetra_INST_CUDA)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    CrsMatrix_WithGraph_Cuda
    SOURCES
      CrsMatrix_WithGraph_Cuda.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    ARGS ${ARGS}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF ()

ASSERT_DEFINED (Tpetra_INST_HIP)
IF (Tpetra_INST_HIP)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    CrsMatrix_WithGraph_HIP
    SOURCES
      CrsMatrix_WithGraph_HIP.cpp
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
      CrsMatrix_WithGraph_OpenMP.cpp
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
      CrsMatrix_WithGraph_Serial.cpp
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
      CrsMatrix_WithGraph_Threads.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    ARGS ${ARGS}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF ()


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_ReplaceDomainMapAndImporter
  SOURCES
    CrsMatrix_ReplaceDomainMapAndImporter.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_ReplaceRangeMapAndExporter
  SOURCES
    CrsMatrix_ReplaceRangeMapAndExporter.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )


TRIBITS_COPY_FILES_TO_BINARY_DIR(CrsMatrixCopyFiles1
  SOURCE_FILES west0067.rua mhd1280b.cua
  EXEDEPS CrsMatrix_UnitTests
  )

# This test only makes sense for > 1 processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_NonlocalSumInto
  SOURCES
    CrsMatrix_NonlocalSumInto.cpp
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
    CrsMatrix_NonlocalSumInto_Ignore.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS
  COMM mpi
  NUM_MPI_PROCS 2-4
  STANDARD_PASS_OUTPUT
  )

# This test only makes sense for exactly 2 MPI processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_Bug5978
  SOURCES
  Bug5978.cpp
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
  Bug6069_1.cpp
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
  Bug6069_2.cpp
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
  Bug6171.cpp
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_Bug8794
  SOURCES
  Bug8794.cpp
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_ReplaceLocalValues
  SOURCES
    CrsMatrix_ReplaceLocalValues.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_ReplaceDiagonal
  SOURCES
    CrsMatrix_ReplaceDiagonal.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_MultipleFillCompletes
  SOURCES
    CrsMatrix_MultipleFillCompletes.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_ReindexColumns
  SOURCES
    CrsMatrix_ReindexColumns.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_TransformValues
  SOURCES
    CrsMatrix_TransformValues.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_GetRowCopy
  SOURCES
    CrsMatrix_GetRowCopy.cpp
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
    assembleElement.cpp
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

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  applyDirichlet
  SOURCES
    applyDirichlet.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_StaticImportExport
  SOURCES
    CrsMatrix_StaticImportExport.cpp
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  sumIntoStaticProfileExtraSpace
  SOURCES sumIntoStaticProfileExtraSpace.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnpackMerge
  SOURCES
    UnpackMerge.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_Bug7745
  SOURCES
  Bug7745.cpp
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS
  COMM mpi serial
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_Bug8447
  SOURCES
  Bug8447.cpp
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_RemoveEmptyProcesses
  SOURCES
    CrsMatrix_removeEmptyProcesses.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2-4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_MatvecH2D
  SOURCES
    CrsMatrix_MatvecH2D.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_MatvecFence
  SOURCES
    CrsMatrix_MatvecFence.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  STANDARD_PASS_OUTPUT
  )


if (
  # supported TPLs
  (
    (Tpetra_INST_CUDA AND TPL_ENABLE_CUSPARSE ) OR
    (Tpetra_INST_HIP AND TPL_ENABLE_ROCSPARSE)
  )

  AND

  # supported type combos
  (
    (Tpetra_INST_DOUBLE OR Tpetra_INST_FLOAT)
  )
)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    CrsMatrix_ApplyUsesTPLs
    SOURCES
      CrsMatrix_ApplyUsesTPLs.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
endif()

SET(TIMING_INSTALLS "")

INSTALL(TARGETS ${TIMING_INSTALLS}
        RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")

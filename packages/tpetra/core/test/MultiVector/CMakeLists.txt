TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVector_UnitTests
  SOURCES
    MultiVector_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  ${MAXNP}
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Bug9856
  SOURCES
    Bug9856.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Bug9583
  SOURCES
    Bug9583.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Bug7758
  SOURCES
    Bug7758.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Bug7745
  SOURCES
    Bug7745.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Bug5474
  SOURCES
    Bug5474.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVector_LocalViewTests
  SOURCES
    MultiVector_LocalViewTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  ${MAXNP}
  )

# This test requires exactly 2 MPI processes.  It could be rewritten
# so that it will run with more, but then the test would just waste
# the remaining processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MV_subViewSomeZeroRows
  SOURCES
    MV_subViewSomeZeroRows.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

# This test only uses one process, and does not require or use MPI,
# but it still works if MPI is enabled.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MV_Kokkos_DualView_subview
  SOURCES
  MV_Kokkos_DualView_subview.cpp
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

# This test comes from zoltan2/test/temp/multivectorTest.cpp.  It
# mimics the recursive bisection algorithm in Zoltan2, and times data
# migration.
#
# The test itself does not refer to MPI.  Hoewver, it only makes sense
# to run this test in an MPI build.  Otherwise, no data migration is
# necessary.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  rcb
  SOURCES rcb.cpp ${MV_RCB_OBJECTS_CUDA}
  COMM mpi
  STANDARD_PASS_OUTPUT
  ${MAXNP}
  )

# This test relates to Github Issue #358.
# It does not need multiple MPI processes,
# but it still works if MPI is enabled.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVector_get2dView
  SOURCES
    MultiVector_get2dView.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

# This test relates to Github Issue #567.
# It does not need multiple MPI processes,
# but it still works if MPI is enabled.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVector_TwoArgRandomize
  SOURCES
    MultiVector_TwoArgRandomize.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

# This test only uses one process, and does not need MPI,
# but it still works if MPI is enabled.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Issue364
  SOURCES
    Issue364.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

# This test only uses one process, and does not need MPI,
# but it should still work if MPI is enabled.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Issue46
  SOURCES
    Issue46.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVector_MicroBenchmark
  SOURCES
    MultiVector_MicroBenchmark.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  ${MAXNP}
  )

ASSERT_DEFINED( ${PACKAGE_NAME}_ENABLE_TeuchosNumerics )
IF (${PACKAGE_NAME}_ENABLE_TeuchosNumerics)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    deep_copy_SerialDenseMatrix_to_MultiVector
    SOURCES
      deep_copy_from_SDM.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    deep_copy_MultiVector_to_SerialDenseMatrix
    SOURCES
      deep_copy_from_SDM.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )
ENDIF ()

# Regression test for Trilinos GitHub Issue #4639.
# It must run with at least 2 MPI processes to exercise
# the fix for that issue.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MV_reduce_strided
  SOURCES
    MV_reduce_strided.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 2-4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Vector_offsetViewCtor
  SOURCES
    Vector_offsetViewCtor.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  aliased_deep_copy
  SOURCES
    aliased_deep_copy.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# The mixed node tests only make sense if there is more than 1 MPI rank
ASSERT_DEFINED (Tpetra_INST_CUDA)
ASSERT_DEFINED (Tpetra_INST_SERIAL)
IF (Tpetra_INST_CUDA AND Tpetra_INST_SERIAL)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MultiVector_MixedNode
    SOURCES
      MultiVector_MixedNode.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    ARGS ${ARGS}
    COMM mpi
    NUM_MPI_PROCS 2-4
    STANDARD_PASS_OUTPUT
    )
ENDIF ()

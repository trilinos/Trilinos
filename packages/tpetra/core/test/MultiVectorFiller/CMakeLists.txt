
ASSERT_DEFINED (Tpetra_ENABLE_Kokkos_Refactor)
IF (Tpetra_ENABLE_Kokkos_Refactor)
  # MultiVectorFiller may use Kokkos data structures.
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MultiVectorFiller
    SOURCES MultiVectorFiller.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF ()


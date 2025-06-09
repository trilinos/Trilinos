TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SingletonFiltering_UnitTests
  SOURCES
    SingletonFiltering_UnitTests.cpp
    SingletonFiltering_TestUtils.hpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SingletonFiltering_np1_UnitTests
  SOURCES
    SingletonFiltering_Parallel_UnitTests.cpp
    SingletonFiltering_TestUtils.hpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SingletonFiltering_np4_UnitTests
  SOURCES
    SingletonFiltering_Parallel_UnitTests.cpp
    SingletonFiltering_TestUtils.hpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(SingletonFilteringFileCopies
  SOURCE_FILES
    SF1_Matrix_Original.mm
    SF1_LHS_Original.mm
    SF1_RHS_Original.mm
    SF1_Matrix_Reduced.mm
    SF1_LHS_Reduced.mm
    SF1_RHS_Reduced.mm
    SF1_Matrix_Reduced_Map.mm
    SF2_Matrix_Original.mm
    SF2_LHS_Original.mm
    SF2_RHS_Original.mm
    SF2_Matrix_Reduced.mm
    SF2_LHS_Reduced.mm
    SF2_RHS_Reduced.mm
  EXEDEPS
    SingletonFiltering_UnitTests
    SingletonFiltering_np1_UnitTests
    SingletonFiltering_np4_UnitTests
  )

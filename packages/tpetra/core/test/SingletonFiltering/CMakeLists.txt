TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SingletonFiltering
  SOURCES
    SingletonFiltering_UnitTests.cpp
    SingletonFiltering_TestUtils.hpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE(
  SingletonFiltering_Parallel
  SOURCES
    SingletonFiltering_Parallel_UnitTests.cpp
    SingletonFiltering_TestUtils.hpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  )

TRIBITS_ADD_TEST(
  SingletonFiltering_Parallel
  NAME SingletonFiltering_Parallel
  ARGS ""
  COMM mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  SingletonFiltering_Parallel
  NAME SingletonFiltering_Parallel
  ARGS ""
  COMM mpi
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
    SF1_Solution_Reduced.mm
    SF1_Solution.mm
    SF2_Matrix_Original.mm
    SF2_LHS_Original.mm
    SF2_RHS_Original.mm
    SF2_Matrix_Reduced.mm
    SF2_LHS_Reduced.mm
    SF2_RHS_Reduced.mm
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SingletonFilteringTransform_UnitTests
  SOURCES
    SingletonFilteringTransform_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(SingletonFilteringTransformFileCopies
  SOURCE_FILES SF1_Matrix_Original.mm SF1_LHS_Original.mm SF1_RHS_Original.mm SF1_Matrix_Reduced.mm SF1_LHS_Reduced.mm SF1_RHS_Reduced.mm SF1_Matrix_Reduced_Map.mm
  EXEDEPS SingletonFilteringTransform_UnitTests
  )

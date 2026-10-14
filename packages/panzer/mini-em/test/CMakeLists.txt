
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TensorConductivity
  SOURCES
     tensor_conductivity.cpp
     ${PANZER_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

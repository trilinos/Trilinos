TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CommandLineProcessor_test  
  SOURCES cxx_main.cpp 
  ARGS --verbose
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CLP_Remove_Quotes
  SOURCES Teuchos_CLP_Remove_Quotes.cpp
  ARGS "--method=TEST"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

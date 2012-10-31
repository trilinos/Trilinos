
INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ConstNonconstObjectContainer_UnitTest  
  SOURCES
    ConstNonconstObjectContainer.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi  
  NUM_MPI_PROCS 1
  )

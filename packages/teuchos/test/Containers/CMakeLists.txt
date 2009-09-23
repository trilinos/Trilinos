
INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ConstNonconstObjectContainer_UnitTest  
  SOURCES
    ConstNonconstObjectContainer.cpp
    ../UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi  
  NUM_MPI_PROCS 1
  )

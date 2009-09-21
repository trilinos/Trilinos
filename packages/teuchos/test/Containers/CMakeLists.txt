
INCLUDE(PackageAddExecutableAndTest)



PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ConstNonconstObjectContainer_UnitTest  
  SOURCES
    ConstNonconstObjectContainer.cpp
    ../UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi  
  NUM_MPI_PROCS 1
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Containers_test  
  SOURCES cxx_main.cpp
  COMM serial mpi  
  XHOSTTYPE Windows
  NUM_MPI_PROCS 1
  )

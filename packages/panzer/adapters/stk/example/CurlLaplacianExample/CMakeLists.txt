INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddAdvancedTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(ASSEMBLY_EXAMPLE_SOURCES
  main.cpp
  )
  
TRIBITS_ADD_EXECUTABLE(
  CurlLaplacianExample 
  SOURCES ${ASSEMBLY_EXAMPLE_SOURCES}
  )

TRIBITS_ADD_ADVANCED_TEST(
  CurlLaplacianExample
  TEST_0 EXEC CurlLaplacianExample
    ARGS --use-epetra
    PASS_REGULAR_EXPRESSION "ALL PASSED: Epetra"
    NUM_MPI_PROCS 4
  TEST_1 EXEC CurlLaplacianExample
    ARGS --use-tpetra
    PASS_REGULAR_EXPRESSION "ALL PASSED: Tpetra"
    NUM_MPI_PROCS 4
  COMM mpi
  )


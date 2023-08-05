TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

IF (${PACKAGE_NAME}_ENABLE_mpi_advance)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MpiAdvance_AllToAllV
    SOURCES
      AllToAllV.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM mpi
    NUM_MPI_PROCS 2
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MpiAdvance_NeighborAllToAllV
    SOURCES
      NeighborAllToAllV.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM mpi
    NUM_MPI_PROCS 2
    STANDARD_PASS_OUTPUT
    )
    
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MpiAdvance_AllToAll
    SOURCES
      AllToAll.cpp
    COMM mpi
    NUM_MPI_PROCS 2
    STANDARD_PASS_OUTPUT
    )
    
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MpiAdvance_NeighborAllToAllW
    SOURCES
      NeighborAllToAllW.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM mpi
    NUM_MPI_PROCS 2
    STANDARD_PASS_OUTPUT
    )
    
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MpiAdvance_NeighborPartLocalityAllToAllV
    SOURCES
      NeighborPartLocalityAllToAllV.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM mpi
    NUM_MPI_PROCS 2
    STANDARD_PASS_OUTPUT
    )
    
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MpiAdvance_NeighborLocalityAllToAllV
    SOURCES
      NeighborLocalityAllToAllV.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM mpi
    NUM_MPI_PROCS 2
    STANDARD_PASS_OUTPUT
    )
ENDIF()
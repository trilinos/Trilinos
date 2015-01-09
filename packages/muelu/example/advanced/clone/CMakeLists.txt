INCLUDE(XpetraAddTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Clone example
#
IF(${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND TpetraClassic_ENABLE_ThrustGPUNode)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Clone-Tpetra
    SOURCES Clone-Tpetra.cpp
    NUM_MPI_PROCS 1
    COMM serial mpi
    )

ENDIF()

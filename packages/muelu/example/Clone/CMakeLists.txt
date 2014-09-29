INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Executables
#

IF(${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND
KokkosClassic_ENABLE_Thrust)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Clone-Tpetra
    SOURCES Clone-Tpetra.cpp
    NUM_MPI_PROCS 1
    COMM serial mpi
    DEPLIBS muelu xpetra muelu-adapters
    )

ENDIF()

INCLUDE(TribitsAddExecutableAndTest)

IF(${PACKAGE_NAME}_ENABLE_Xpetra)
    
IF (NOT Xpetra_Epetra_NO_32BIT_GLOBAL_INDICES)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      XpetraMap
      SOURCES XpetraMap.cpp
      NUM_MPI_PROCS 4
      COMM serial mpi
      )
ENDIF()


IF (NOT Xpetra_Epetra_NO_64BIT_GLOBAL_INDICES)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      XpetraMap_LL
      SOURCES XpetraMap.cpp
      NUM_MPI_PROCS 4
      COMM serial mpi
      DEFINES -DXPETRA_TEST_USE_LONGLONG_GO
      )
ENDIF()

ENDIF()

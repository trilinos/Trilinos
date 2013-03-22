INCLUDE(TribitsAddExecutableAndTest)

IF(${PACKAGE_NAME}_ENABLE_Xpetra)
    
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      XpetraMap
      SOURCES XpetraMap.cpp
      NUM_MPI_PROCS 4
      COMM serial mpi
      )

ENDIF()

ADD_SUBDIRECTORY(VerySimple)

INCLUDE(TribitsAddExecutableAndTest)

IF(${PACKAGE_NAME}_ENABLE_Xpetra AND NOT Trilinos_ENABLE_EXPLICIT_INSTANTIATION))
    
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      bug_5851
      SOURCES bug_5851.cpp
      NUM_MPI_PROCS 4
      COMM serial mpi
      )

ENDIF()


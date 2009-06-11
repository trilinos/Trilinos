
INCLUDE(PackageAddExecutableAndTest)

IF(${PACKAGE_NAME}_ENABLE_Galeri)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Test_Epetra_RowMatrix
    SOURCES cxx_main.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
    )
ENDIF()

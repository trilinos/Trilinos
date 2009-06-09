INCLUDE(PackageAddExecutableAndTest)

IF(Trilinos_ENABLE_Epetra AND Trilinos_ENABLE_EpetraExt AND
   Trilinos_ENABLE_Pamgen)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Example_01
    SOURCES example_01.cpp
    DEPLIBS epetra epetraext pamgen
    ARGS "10 10 10 false 1.0 10.0 0.0 1.0 -1.0 1.0 -1.0 1.0"
    NUM_MPI_PROCS 1
    ADD_DIR_TO_NAME
    )
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Example_02
    SOURCES example_02.cpp
    DEPLIBS epetra epetraext pamgen
    ARGS "10 10 10 false 1.0 10.0 0.0 1.0 -1.0 1.0 -1.0 1.0"
    NUM_MPI_PROCS 1
    ADD_DIR_TO_NAME
    )
ENDIF()

IF(Trilinos_ENABLE_Epetra AND Trilinos_ENABLE_EpetraExt AND
   Trilinos_ENABLE_Pamgen AND Trilinos_ENABLE_AztecOO AND
    Trilinos_ENABLE_Amesos AND Trilinos_ENABLE_ML)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Example_Poisson
    SOURCES example_Poisson.cpp
    DEPLIBS epetra epetraext pamgen aztecoo amesos ml
    ARGS "10 10 10"
    NUM_MPI_PROCS 1
    ADD_DIR_TO_NAME
    )
ENDIF()


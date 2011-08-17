INCLUDE(PackageAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#######################################################################################
#######################################################################################

### MatrixVectorCheckerExample:

# TODO: IF (Xpetra_ENABLE_Epetra) or IF (Xpetra_ENABLE_Tpetra)

PACKAGE_ADD_EXECUTABLE(
  MatrixVectorCheckerExample
  SOURCES MatrixVectorCheckerExample.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra
  )

IF (Xpetra_ENABLE_Epetra)
  PACKAGE_ADD_TEST(
    MatrixVectorCheckerExample
    NAME "MatrixVectorChecker Epetra"
    ARGS "--linAlgebra=0"
    COMM serial mpi
    # NUM_MPI_PROCS 2
    # PASS_REGULAR_EXPRESSION "Example Passed!"
    )
ENDIF()

IF (Xpetra_ENABLE_Tpetra)
  PACKAGE_ADD_TEST(
    MatrixVectorCheckerExample
    NAME "MatrixVectorChecker Epetra"
    ARGS "--linAlgebra=1"
    COMM serial mpi
    # NUM_MPI_PROCS 2
    # PASS_REGULAR_EXPRESSION "Example Passed!"
    )
ENDIF()

#######################################################################################
#######################################################################################

### read_MatrixMarket:

# Note: this example don't use xpetra. Only need Epetra + Galeri.
IF (${PACKAGE_NAME}_ENABLE_Epetra) 
  IF (${PACKAGE_NAME}_ENABLE_Galeri)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      read_MatrixMarket
      SOURCES MueLu_read_MatrixMarket.cpp
      COMM serial
      DEPLIBS muelu muelu-xpetra
      )
  ENDIF()
ENDIF()

#PACKAGE_COPY_FILES_TO_BINARY_DIR(read_MatrixMarket_example_copyfiles
#  SOURCE_FILES rgen.m
#)

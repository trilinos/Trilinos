INCLUDE(TribitsAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#######################################################################################
#######################################################################################

### MatrixVectorCheckerExample:

# TODO: IF (Xpetra_ENABLE_Epetra) or IF (Xpetra_ENABLE_Tpetra)
# FIXME: does not work in parallel
IF (${PACKAGE_NAME}_ENABLE_Galeri)
  TRIBITS_ADD_EXECUTABLE(
    MatrixVectorCheckerExample
    SOURCES MatrixVectorCheckerExample.cpp
    COMM serial
    DEPLIBS muelu xpetra
    )

  IF (${PACKAGE_NAME}_ENABLE_Epetra)
    TRIBITS_ADD_TEST(
      MatrixVectorCheckerExample
      #NAME "MatrixVectorChecker Epetra"
      NAME "MatrixVectorChecker"
      ARGS "--linAlgebra=0"
      COMM serial
      # NUM_MPI_PROCS 2
      # PASS_REGULAR_EXPRESSION "Example Passed!"
      )
  ENDIF()
  
  IF (${PACKAGE_NAME}_ENABLE_Tpetra)
    TRIBITS_ADD_TEST(
      MatrixVectorCheckerExample
      #NAME "MatrixVectorChecker Epetra"
      NAME "MatrixVectorChecker"
      ARGS "--linAlgebra=1"
      COMM serial
      # NUM_MPI_PROCS 2
      # PASS_REGULAR_EXPRESSION "Example Passed!"
      )
  ENDIF()

ENDIF()

#######################################################################################
#######################################################################################

### read_MatrixMarket:

# Note: this example don't use xpetra. Only need Epetra + Galeri.
IF (${PACKAGE_NAME}_ENABLE_Epetra) 
  IF (${PACKAGE_NAME}_ENABLE_Galeri)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      read_MatrixMarket
      SOURCES MueLu_read_MatrixMarket.cpp
      COMM serial mpi
      DEPLIBS muelu xpetra
      )
  ENDIF()
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(read_MatrixMarket_example_copyfiles
  SOURCE_FILES mat_example.mm #rgen.m
)

INCLUDE(PackageAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#######################################################################################
#######################################################################################

# MatrixDriver

# TODO: IF (Cthulhu_ENABLE_Epetra) or IF (Cthulhu_ENABLE_Tpetra)
PACKAGE_ADD_EXECUTABLE(
  Cthulhu_MatrixDriver
  SOURCES MueLu_MatrixDriver.cpp
  COMM serial mpi
  DEPLIBS muelu-cthulhu
  )

IF (Cthulhu_ENABLE_Epetra)
  PACKAGE_ADD_TEST(
    Cthulhu_MatrixDriver
    NAME "MatrixDriver Epetra"
    ARGS "--linAlgebra=0"
    COMM serial mpi
    # NUM_MPI_PROCS 2
    # PASS_REGULAR_EXPRESSION "Example Passed!"
    )
ENDIF()

IF (Cthulhu_ENABLE_Tpetra)
  PACKAGE_ADD_TEST(
    Cthulhu_MatrixDriver
    NAME "MatrixDriver Tpetra"
    ARGS "--linAlgebra=1"
    COMM serial mpi
    # NUM_MPI_PROCS 2
    # PASS_REGULAR_EXPRESSION "Example Passed!"
    )
ENDIF()


#######################################################################################
#######################################################################################

# TODO: IF (Cthulhu_ENABLE_Epetra) or IF (Cthulhu_ENABLE_Tpetra)
PACKAGE_ADD_EXECUTABLE(
  Cthulhu_MatrixDriverTests
  SOURCES MueLu_MatrixDriverTests.cpp
  COMM serial mpi
  DEPLIBS muelu-cthulhu
  )

IF (Cthulhu_ENABLE_Epetra)
  PACKAGE_ADD_TEST(
    Cthulhu_MatrixDriverTests
    NAME "MatrixDriverTests Epetra"
    ARGS "--linAlgebra=0"
    COMM serial mpi
    # NUM_MPI_PROCS 2
    # PASS_REGULAR_EXPRESSION "Example Passed!"
    )
ENDIF()

IF (Cthulhu_ENABLE_Tpetra)
  PACKAGE_ADD_TEST(
    Cthulhu_MatrixDriverTests
    NAME "MatrixDriverTests Tpetra"
    ARGS "--linAlgebra=1"
    COMM serial mpi
    # NUM_MPI_PROCS 2
    # PASS_REGULAR_EXPRESSION "Example Passed!"
    )
ENDIF()



INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../utils)

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    NOX_Operators
    SOURCES
    test.C 
    Laplace2D.C
    Laplace2D.H
    DEPLIBS noxtestutils
    COMM serial mpi
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Test passed!"
    )
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    BroydenOp
    SOURCES
    test_BroydenOp.C
    DEPLIBS noxtestutils
    COMM serial mpi
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Test passed!"
    )
  
ENDIF()

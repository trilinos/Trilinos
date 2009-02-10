
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_EpetraExt AND NOX_ENABLE_LOCA AND TPL_ENABLE_MPI)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_MultiPointTcubed
    SOURCES 
    Problem_Interface_MP.H
    Problem_Interface_MP.C
    multipointTcubed.C
    COMM mpi
    DEPLIBS noxtestutils locaepetratestproblems
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

ENDIF()

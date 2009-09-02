
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_LOCA)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_Brusselator
    SOURCES 
    Basis.H
    Basis.C
    Brusselator.H
    Brusselator.C
    Problem_Interface.H
    Problem_Interface.C
    brussCont.C
    COMM serial mpi
    DEPLIBS noxtestutils
    NUM_MPI_PROCS 2
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  IF(NOX_ENABLE_EpetraExt)

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      LOCA_BrusselatorHopf
      SOURCES 
      Basis.H
      Basis.C
      Brusselator.H
      Brusselator.C
      Problem_Interface.H
      Problem_Interface.C
      brussHopfCont.C
      COMM serial mpi
      NUM_MPI_PROCS 2
      DEPLIBS noxtestutils
      ARGS -v
      PASS_REGULAR_EXPRESSION "All tests passed"
      )

  ENDIF()

ENDIF()

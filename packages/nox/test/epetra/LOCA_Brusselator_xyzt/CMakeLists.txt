

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../../utils)

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_EpetraExt AND NOX_ENABLE_LOCA)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_BrussXYZT_Ifpack
    SOURCES
    Basis.H
    Basis.C
    Brusselator.H
    Brusselator.C
    Problem_Interface.H
    Problem_Interface.C
    brussXYZT_ifpack.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    DEPLIBS noxtestutils
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_BrussXYZT_Global
    SOURCES
    Basis.H
    Basis.C
    Brusselator.H
    Brusselator.C
    Problem_Interface.H
    Problem_Interface.C
    brussXYZT_global.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    DEPLIBS noxtestutils
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_BrussXYZT_Sequential
    SOURCES
    Basis.H
    Basis.C
    Brusselator.H
    Brusselator.C
    Problem_Interface.H
    Problem_Interface.C
    brussXYZT_sequential.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    DEPLIBS noxtestutils
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_BrussXYZT_SequentialIPS
    SOURCES
    Basis.H
    Basis.C
    Brusselator.H
    Brusselator.C
    Problem_Interface.H
    Problem_Interface.C
    brussXYZT_sequential_IPS.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    DEPLIBS noxtestutils
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_BrussXYZT_SequentialOPS
    SOURCES
    Basis.H
    Basis.C
    Brusselator.H
    Brusselator.C
    Problem_Interface.H
    Problem_Interface.C
    brussXYZT_sequential_OPS.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    DEPLIBS noxtestutils
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_BrussXYZT_Parallel
    SOURCES
    Basis.H
    Basis.C
    Brusselator.H
    Brusselator.C
    Problem_Interface.H
    Problem_Interface.C
    brussXYZT_parallel.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    DEPLIBS noxtestutils
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  # Commented out from original Makefile.am
  # breaks mpirun
#  TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    LOCA_BrussXYZT_Parareal
#    SOURCES
#    Basis.H
#    Basis.C
#    Brusselator.H
#    Brusselator.C
#    Problem_Interface.H
#    Problem_Interface.C
#    brussXYZT_parareal.C
#    COMM serial mpi
#    DEPLIBS noxtestutils
#    ARGS -v
#    PASS_REGULAR_EXPRESSION "All tests passed"
#    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_BrussXYZT_BlockDiagonal
    SOURCES
    Basis.H
    Basis.C
    Brusselator.H
    Brusselator.C
    Problem_Interface.H
    Problem_Interface.C
    brussXYZT_blockdiagonal.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    DEPLIBS noxtestutils
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_BrussXYZT_bdsdt
    SOURCES
    Basis.H
    Basis.C
    Brusselator.H
    Brusselator.C
    Problem_Interface.H
    Problem_Interface.C
    brussXYZT_bdsdt.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    DEPLIBS noxtestutils
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

ENDIF()



TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_ModalSolverUtils
  SOURCES cxx_main.cpp
  ARGS
  COMM serial mpi
  )

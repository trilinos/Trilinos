TRIBITS_ADD_EXECUTABLE(
  Tpetra_TD_Gen_example
  SOURCES TraceMinDavidsonGeneralizedEx.cpp
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE(
  Tpetra_TD_Fiedler_example
  SOURCES TraceMinDavidsonLaplacianEx.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  Tpetra_TD_Trans_example
  SOURCES TraceMinDavidsonSpecTransEx.cpp
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE(
  Tpetra_TD_UserOp_example
  SOURCES TraceMinDavidsonUserOpEx.cpp
  COMM serial mpi
)

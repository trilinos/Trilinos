TRIBITS_ADD_EXECUTABLE(
  Tpetra_LOBPCG_CustomStatusTest
  SOURCES LOBPCGCustomStatusTest.cpp
  COMM serial mpi
)

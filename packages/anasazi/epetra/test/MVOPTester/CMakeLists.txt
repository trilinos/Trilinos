
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Epetra_MVOPTester
  SOURCES cxx_main.cpp
  COMM serial mpi
  ARGS "-v"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE(
  Epetra_SpecializedMVOPTester
  SOURCES cxx_specialized_main.cpp
  COMM serial mpi
  )

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  Epetra_SpecializedMVOPTester
#  SOURCES cxx_specialized_main.cpp
#  COMM serial mpi
#  ARGS "-v"
#  STANDARD_PASS_OUTPUT
#  )

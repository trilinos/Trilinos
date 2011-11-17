INCLUDE(TribitsAddExecutableAndTest)

IF (Anasazi_ENABLE_Epetra)
TRIBITS_ADD_EXECUTABLE(
  Thyra_MVOPTester_example
  SOURCES MVOPTesterEx.cpp
  COMM serial mpi
)
ENDIF()

INCLUDE(PackageAddExecutableAndTest)

IF (Anasazi_ENABLE_Epetra)
PACKAGE_ADD_EXECUTABLE(
  Thyra_MVOPTester_example
  SOURCES MVOPTesterEx.cpp
  COMM serial mpi
)
ENDIF()

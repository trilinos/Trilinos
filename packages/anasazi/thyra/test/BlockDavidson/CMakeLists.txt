
INCLUDE(PackageAddExecutableAndTest)

IF(Anasazi_ENABLE_Epetra)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    BlockDavidsonThyra_test
    SOURCES cxx_main.cpp
    ARGS 
    COMM serial mpi
    )
ENDIF()




IF (NOT Xpetra_Epetra_NO_32BIT_GLOBAL_INDICES)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Cloner_UnitTests
    SOURCES
    Cloner_UnitTests
      ../Xpetra_UnitTests
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF()

ASSERT_DEFINED(Belos_ENABLE_EpetraExt)
ASSERT_DEFINED(Belos_ENABLE_Triutils)
ASSERT_DEFINED(Belos_ENABLE_ML)
IF (Belos_ENABLE_Triutils AND Belos_ENABLE_EpetraExt AND Belos_ENABLE_ML)
  TRIBITS_ADD_EXECUTABLE(
    pcpg_epetraex
    SOURCES test_pcpg_epetraex.cpp
    COMM serial mpi
    )
ENDIF()

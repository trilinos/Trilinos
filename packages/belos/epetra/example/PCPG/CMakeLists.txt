ASSERT_DEFINED(Belos_ENABLE_EpetraExt)
ASSERT_DEFINED(Belos_ENABLE_Triutils)
ASSERT_DEFINED(Belos_ENABLE_ML)
IF (Belos_ENABLE_Triutils AND Belos_ENABLE_EpetraExt AND Belos_ENABLE_ML)
  TRIBITS_ADD_EXECUTABLE(
    PCPG_Epetra_File_Ex
    SOURCES PCPGEpetraExFile.cpp
    COMM serial mpi
    )
ENDIF()

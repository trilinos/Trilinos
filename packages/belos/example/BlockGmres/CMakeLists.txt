
INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Copy_Files_To_Binary_Dir)


ASSERT_DEFINED(Trilinos_ENABLE_Triutils)
ASSERT_DEFINED(Trilinos_ENABLE_EpetraExt)
ASSERT_DEFINED(Trilinos_ENABLE_Ifpack)
IF(Trilinos_ENABLE_Triutils AND Trilinos_ENABLE_EpetraExt AND Trilinos_ENABLE_Ifpack)

  
  TRILINOS_ADD_EXECUTABLE_AND_TEST(
    BlockPrecGmresEpetraExFile
    COMM serial mpi
    SOURCES BlockPrecGmresEpetraExFile.cpp
    )

  COPY_FILES_TO_BINARY_DIR(CopyExampleBlockGmresBelosFiles
    SOURCE_FILES orsirr1.hb fidap036.hb
    EXEDEPS BlockPrecGmresEpetraExFile
    )


ENDIF()

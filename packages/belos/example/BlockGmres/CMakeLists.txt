
INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Copy_Files_To_Binary_Dir)


IF(Trilinos_ENABLE_Triutils AND Trilinos_ENABLE_EpetraExt 
    AND Trilinos_ENABLE_Ifpack)

  
  TRILINOS_ADD_EXECUTABLE(
    BlockPrecGmresEpetraExFile
    COMM serial mpi
    SOURCES BlockPrecGmresEpetraExFile.cpp
    )

  COPY_FILES_TO_BINARY_DIR(CopyExampleBlockGmresBelosFiles
    SOURCE_FILES orsirr1.hb fidap036.hb
    )


ENDIF()

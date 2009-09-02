
INCLUDE(PackageAddExecutable)
INCLUDE(PackageCopyFilesToBinaryDir)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
IF (${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_EpetraExt)
    
  PACKAGE_ADD_EXECUTABLE(
    BlockGmres_Epetra_File_Ex
    SOURCES BlockGmresEpetraExFile.cpp
    COMM serial mpi
    )

  PACKAGE_ADD_EXECUTABLE(
    PseudoBlockGmres_Epetra_File_Ex
    SOURCES PseudoBlockGmresEpetraExFile.cpp
    COMM serial mpi
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyExampleBlockGmresFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/epetra/example/BlockGmres
    SOURCE_FILES orsirr1.hb
    EXEDEPS BlockGmres_Epetra_File_Ex
    )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF(${PACKAGE_NAME}_ENABLE_Ifpack)

    PACKAGE_ADD_EXECUTABLE(
      BlockPrecGmres_Epetra_File_Ex
      SOURCES BlockPrecGmresEpetraExFile.cpp
      COMM serial mpi
      )

    PACKAGE_ADD_EXECUTABLE(
      PseudoBlockPrecGmres_Epetra_File_Ex
      SOURCES PseudoBlockPrecGmresEpetraExFile.cpp
      COMM serial mpi
      )

    PACKAGE_ADD_EXECUTABLE(
      BlockFlexPrecGmres_Epetra_File_Ex
      SOURCES BlockFlexGmresEpetraExFile.cpp
      COMM serial mpi
      )

  ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_EpetraExt)

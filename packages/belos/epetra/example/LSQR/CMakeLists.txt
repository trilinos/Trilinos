
INCLUDE(PackageAddExecutable)
INCLUDE(PackageCopyFilesToBinaryDir)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
IF (${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_EpetraExt)
    
  PACKAGE_ADD_EXECUTABLE(
    LSQR_Epetra_File_Ex
    SOURCES LSQREpetraExFile.cpp
    COMM serial mpi
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyExampleLSQRFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/epetra/example/LSQR
    SOURCE_FILES orsirr1_scaled.hb
    EXEDEPS LSQR_Epetra_File_Ex
    )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF(${PACKAGE_NAME}_ENABLE_Ifpack)

    PACKAGE_ADD_EXECUTABLE(
      PrecLSQR_Epetra_File_Ex
      SOURCES PrecLSQREpetraExFile.cpp
      COMM serial mpi
      )



  ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_EpetraExt)

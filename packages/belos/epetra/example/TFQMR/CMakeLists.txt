
INCLUDE(PackageAddExecutable)
INCLUDE(PackageCopyFilesToBinaryDir)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF (${PACKAGE_NAME}_ENABLE_Ifpack)

  PACKAGE_ADD_EXECUTABLE(
    tfqmr_ex
    SOURCES TFQMREpetraExFile.cpp createEpetraProblem.cpp
    COMM serial mpi
  )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyExampleTFQMRFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/epetra/example/BlockGmres
    SOURCE_FILES orsirr1.hb 
    EXEDEPS tfqmr_ex
    )

  ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)

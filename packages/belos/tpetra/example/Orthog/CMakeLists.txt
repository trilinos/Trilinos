
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Tpetra)
IF(${PACKAGE_NAME}_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE(
    TpetraOrthogEx
    SOURCES TpetraOrthogEx.cpp
    COMM serial mpi
    )

  ASSERT_DEFINED(Belos_SOURCE_DIR)
  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_OrthogEx_CopyFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/tpetra/example/Orthog
    SOURCE_FILES orsirr_1.mtx
    EXEDEPS TpetraOrthogEx
    )
ENDIF()

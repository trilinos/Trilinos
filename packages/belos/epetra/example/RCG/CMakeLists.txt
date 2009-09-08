
INCLUDE(PackageAddExecutable)
INCLUDE(PackageCopyFilesToBinaryDir)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_EpetraExt)
    
  PACKAGE_ADD_EXECUTABLE(
    RCG_Epetra_Ex_File
    SOURCES RCGEpetraExFile.cpp createEpetraProblem.cpp
    COMM serial mpi
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyExampleRCGFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/epetra/example/RCG
    SOURCE_FILES bcsstk14.hb
    EXEDEPS RCG_Epetra_Ex_File
    )

ENDIF()

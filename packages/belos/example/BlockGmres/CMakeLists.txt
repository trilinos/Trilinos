
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)


ASSERT_DEFINED(Belos_ENABLE_Triutils)
ASSERT_DEFINED(Belos_ENABLE_EpetraExt)
ASSERT_DEFINED(Belos_ENABLE_Ifpack)
IF(Belos_ENABLE_Triutils AND Belos_ENABLE_EpetraExt AND Belos_ENABLE_Ifpack)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    BlockPrecGmresEpetraExFile
    SOURCES BlockPrecGmresEpetraExFile.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyExampleBlockGmresBelosFiles
    SOURCE_FILES orsirr1.hb fidap036.hb
    EXEDEPS BlockPrecGmresEpetraExFile
    )

ENDIF()

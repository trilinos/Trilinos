
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  read_test
  SOURCES TestTriutilReads.cpp
  COMM serial mpi
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(ReadTestCopyFiles
  DEST_FILES
    bcsstk01.rsa bcsstk01.mtx bcsstk01.triS Diagonal.rua Diagonal.mtx
    Diagonal.TimD fs_183_4.rua fs_183_4.mtx fs_183_4.triU
    impcol_a.rua impcol_a.mtx impcol_a.triU
  SOURCE_PREFIX "In_"
  EXEDEPS read_test
  )

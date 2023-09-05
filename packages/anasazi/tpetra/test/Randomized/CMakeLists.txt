
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_Randomized_Experimental
  SOURCES cxx_experimental.cpp 
  ARGS 
    "--verbose"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_Randomized_DoubleFloat_Test
  SOURCES cxx_main_doublefloat.cpp
  ARGS 
    "--verbose"
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_Randomized_DoubleFloatCopyFiles
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/testmatrices
  SOURCE_FILES simple.mtx
  EXEDEPS Tpetra_Randomized_Experimental Tpetra_Randomized_DoubleFloat_Test
  )


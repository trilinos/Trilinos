
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_LOBPCG_test
  SOURCES cxx_main.cpp
  ARGS
    "--verbose"
    "--verbose --ortho=ICGS"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_LOBPCG_auxtest
  SOURCES cxx_main_auxtest.cpp
  ARGS
    "--verbose"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_LOBPCG_solvertest
  SOURCES cxx_main_solvertest.cpp
  ARGS
    "--verbose"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_LOBPCG_simpletest
  SOURCES cxx_simple.cpp
  ARGS
    "--verbose"
  COMM serial mpi
  )

IF (Tpetra_INST_COMPLEX_DOUBLE)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_LOBPCG_Complex_test
    SOURCES cxx_main_complex.cpp 
    ARGS 
      "--verbose"
      "--verbose --ortho=ICGS"
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_LOBPCG_ComplexCopyFiles
    SOURCE_DIR ${PACKAGE_SOURCE_DIR}/testmatrices
    SOURCE_FILES mhd1280b.cua mhd1280b.mtx
    EXEDEPS Tpetra_LOBPCG_Complex_test
    )

ENDIF()


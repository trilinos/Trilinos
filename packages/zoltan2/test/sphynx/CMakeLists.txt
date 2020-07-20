TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Sphynx
  SOURCES Test_Sphynx.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Sphynx
  NAME Sphynx_VWeights
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --vertexWeights=3"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Sphynx
  NAME Sphynx_OneProc_VWeights
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --no-distribute --vertexWeights=1"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Sphynx
  NAME Sphynx_NormalizedLaplacian
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --normalized"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Sphynx
  NAME Sphynx_GeneralizedLaplacian
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --generalized"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )


IF(Trilinos_ENABLE_Galeri)
  TRIBITS_ADD_TEST(
    Sphynx
    NAME Sphynx_Galeri
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

  TRIBITS_ADD_TEST(
    Sphynx
    NAME Sphynx_Galeri_VWeights
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
     "--vertexWeights=3"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

ENDIF()

IF(Trilinos_ENABLE_MueLu)
  TRIBITS_ADD_TEST(
    Sphynx
    NAME Sphynx_PolyPrec
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
     "--inputFile=simple --poly"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()


# Command that copies files to the executable directory.
TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_files_for_sphynx_tests
  SOURCE_FILES
    simple.mtx
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../core/data
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_levelsolver
  SOURCES LevelSolver_test.cpp TestUtils.cpp
  COMM serial mpi
  ARGS
    "M01.mtx 4 10"
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi tpetra
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(TestLevelSolverCopyFiles
  SOURCE_FILES M01.mtx M02.mtx P01.p P02.p
  EXEDEPS test_levelsolver
  )

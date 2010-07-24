INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  SimpleSolveExample
  SOURCES SimpleSolve
  ARGS
  COMM serial mpi
)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  SimpleSolveFileExample
  SOURCES SimpleSolve_File
  ARGS
  COMM serial mpi
)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MultipleSolvesFileExample
  SOURCES MultipleSolves_File
  ARGS
  COMM serial mpi
)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  SimpleSolveWithParameters
  SOURCES SimpleSolve_WithParameters
  ARGS
  COMM serial mpi
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(Amesos2_SimpleSolve_File
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/epetra/example/BlockCG
  SOURCE_FILES bcsstk14.hb
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/amesos2/example
  SOURCE_FILES amesos_test_mat.hb
  EXEDEPS SimpleSolveFileExample
  )

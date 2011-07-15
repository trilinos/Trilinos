INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

# PACKAGE_ADD_EXECUTABLE_AND_TEST(
#   SimpleSolveExample
#   SOURCES SimpleSolve
#   ARGS
#   COMM serial mpi
# )

# PACKAGE_ADD_EXECUTABLE_AND_TEST(
#   SimpleSolveFileExample
#   SOURCES SimpleSolve_File
#   ARGS
#   COMM serial mpi
# )

# PACKAGE_ADD_EXECUTABLE_AND_TEST(
#   MultipleSolvesFileExample
#   SOURCES MultipleSolves_File
#   ARGS
#   COMM serial mpi
# )

# PACKAGE_ADD_EXECUTABLE_AND_TEST(
#   SimpleSolveWithParameters
#   SOURCES SimpleSolve_WithParameters
#   ARGS
#   COMM serial mpi
# )

# PACKAGE_ADD_EXECUTABLE_AND_TEST(
#   ComplexSolveWithTranspose
#   SOURCES ComplexSolve_WithTranspose
#   ARGS
#   COMM serial mpi
# )

# PACKAGE_ADD_EXECUTABLE(
#   TwoPartSolveExample
#   SOURCES TwoPartSolve
#   )

PACKAGE_ADD_EXECUTABLE(
  QuickSolve
  SOURCES quick_solve
  )

IF (${PACKAGE_NAME}_ENABLE_EpetraExt)
  PACKAGE_ADD_EXECUTABLE(
    QuickSolveEpetra
    SOURCES quick_solve_epetra
  )
ENDIF()

PACKAGE_COPY_FILES_TO_BINARY_DIR(Amesos2_SimpleSolve_File
  SOURCE_FILES arc130.mtx
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../test/matrices/
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )

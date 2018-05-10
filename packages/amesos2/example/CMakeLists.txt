
IF (${PACKAGE_NAME}_ENABLE_KLU2 AND NOT ${PROJECT_NAME}_ENABLE_DEBUG AND NOT KOKKOS_ENABLE_DEBUG)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SimpleSolveExampleNonContigMap
  SOURCES SimpleSolveNonContigMap
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 2
)
ENDIF()


TRIBITS_ADD_EXECUTABLE(
  SimpleSolveExample
  SOURCES SimpleSolve
  COMM serial mpi
  )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   MultipleSolvesFileExample
#   SOURCES MultipleSolves_File
#   ARGS
#   COMM serial mpi
# )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   SimpleSolveWithParameters
#   SOURCES SimpleSolve_WithParameters
#   ARGS
#   COMM serial mpi
# )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   ComplexSolveWithTranspose
#   SOURCES ComplexSolve_WithTranspose
#   ARGS
#   COMM serial mpi
# )

# TRIBITS_ADD_EXECUTABLE(
#   TwoPartSolveExample
#   SOURCES TwoPartSolve
#   )


IF (Tpetra_INST_INT_INT)
TRIBITS_ADD_EXECUTABLE(
  QuickSolve
  SOURCES quick_solve
  )
ENDIF ()

IF (Tpetra_INST_INT_INT)
 TRIBITS_ADD_EXECUTABLE(
   SimpleSolveFileExample
  SOURCES SimpleSolve_File
 )
ENDIF ()


IF (Tpetra_INST_INT_INT)
IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt)
  TRIBITS_ADD_EXECUTABLE(
    QuickSolveEpetra
    SOURCES quick_solve_epetra
  )
ENDIF()
ENDIF()

IF (Tpetra_INST_INT_INT)
TRIBITS_COPY_FILES_TO_BINARY_DIR(Amesos2_SimpleSolve_File
  SOURCE_FILES arc130.mtx
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../test/matrices/
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )
ENDIF()

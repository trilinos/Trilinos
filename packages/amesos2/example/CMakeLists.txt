
IF (${PACKAGE_NAME}_ENABLE_KLU2 AND NOT ${PROJECT_NAME}_ENABLE_DEBUG AND NOT Kokkos_ENABLE_DEBUG)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SimpleSolveExampleNonContigMap
  SOURCES SimpleSolveNonContigMap.cpp
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 2
)
ENDIF()

IF (Tpetra_INST_DOUBLE)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  GappedMtxGIDs-1proc
  SOURCES GappedMtxGIDs-1proc.cpp
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1
)
ENDIF ()

TRIBITS_COPY_FILES_TO_BINARY_DIR(GappedMtxGIDs_File
  SOURCE_FILES gap-ids-1proc-rowmap-tpetra.mm gap-ids-1procA.mm gap-ids-1procRhs.mm
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../test/matrices/
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )

TRIBITS_ADD_EXECUTABLE(
  SimpleSolveExample
  SOURCES SimpleSolve.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  SimpleSolveTransposeExample
  SOURCES SimpleSolveTranspose.cpp
  COMM serial mpi
  )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   MultipleSolvesFileExample
#   SOURCES MultipleSolves_File.cpp
#   ARGS
#   COMM serial mpi
# )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   SimpleSolveWithParameters
#   SOURCES SimpleSolve_WithParameters.cpp
#   ARGS
#   COMM serial mpi
# )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   ComplexSolveWithTranspose
#   SOURCES ComplexSolve_WithTranspose.cpp
#   ARGS
#   COMM serial mpi
# )

# TRIBITS_ADD_EXECUTABLE(
#   TwoPartSolveExample
#   SOURCES TwoPartSolve.cpp
#   )


TRIBITS_ADD_EXECUTABLE(
  QuickSolve
  SOURCES quick_solve.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  SimpleSolveFileExample
  SOURCES SimpleSolve_File.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  SimpleSolveIRExample
  SOURCES SimpleSolve_IR.cpp
  )

IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt)
  TRIBITS_ADD_EXECUTABLE(
    QuickSolveEpetra
    SOURCES quick_solve_epetra.cpp
  )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(Amesos2_SimpleSolve_File
  SOURCE_FILES arc130.mtx
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../test/matrices/
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )

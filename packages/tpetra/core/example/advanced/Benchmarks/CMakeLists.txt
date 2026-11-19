IF(Tpetra_INST_INT_INT)
  # This example requires GlobalOrdinal = int.
  TRIBITS_ADD_EXECUTABLE(
    import
    SOURCES import.cpp
    COMM mpi
  )
ENDIF()

TRIBITS_ADD_EXECUTABLE(
  vector
  SOURCES vector.cpp
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE(
  blockCrsMatrixMatVec
  SOURCES blockCrsMatrixMatVec.cpp
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE(
  CrsMatrixMatVec
  SOURCES CrsMatrixMatVec.cpp
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE(
  localView
  SOURCES localView.cpp
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE(
  CrsMatrix_sumIntoLocalValues
  SOURCES CrsMatrix_sumIntoLocalValues.cpp
  COMM serial mpi
)

# This benchmark needs at least 2 MPI processes.
# It builds fine without MPI, but there's no point
# in building it for that case.
TRIBITS_ADD_EXECUTABLE(
  CrsMatrixDenseRowUnpack
  SOURCES CrsMatrixDenseRowUnpack.cpp
  COMM mpi
)

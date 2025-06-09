
# for TpetraTestingUtilities.hpp
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

ADD_SUBDIRECTORIES(
  BasicPerfTest
  Behavior
  Blas
  Block
  BugTests
  Comm
  Core
  CrsGraph
  CrsGraphTransposer
  CrsMatrix
  Directory
  Distributor
  FECrsGraph
  FECrsMatrix
  FEMultiVector
  Filter
  HashTable
  ImportExport
  ImportExport2
  inout
  LinearProblem
  Map
  MatrixMatrix
  Merge
  MPIAdvance
  MultiVector
  Node
  Operator
  PerformanceCGSolve
  Sort
  Spaces
  Utils
  RowMatrixTransposer
  Tsqr
  ReindexTransform
  SolverMapTransform
  SingletonFiltering
  )

IF (Tpetra_ENABLE_DEPRECATED_CODE AND ${PACKAGE_NAME}_ENABLE_Epetra)
  ADD_SUBDIRECTORY(EpetraRowMatrix)
ENDIF()

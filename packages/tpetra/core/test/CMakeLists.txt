
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
  CrsMatrix
  Directory
  Distributor
  FEMultiVector
  HashTable
  ImportExport
  ImportExport2
  inout
  Map
  MatrixMatrix
  Merge
  MultiVector
  MultiVectorFiller
  Node
  PerformanceCGSolve
  Sort
  Utils
  RowMatrixTransposer
  TypeStack
  )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Epetra)
IF (${PACKAGE_NAME}_ENABLE_Epetra)
  ADD_SUBDIRECTORY(EpetraRowMatrix)
ENDIF()

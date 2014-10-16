
# for TpetraTestingUtilities.hpp
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

ADD_SUBDIRECTORIES(
  Platform
  Map
  Distributor
  Directory
  ImportExport
  ImportExport2
  MultiVector
  MultiVectorFiller
  BasicPerfTest
  Utils
  CrsGraph
  CrsMatrix
  BlockCrsGraph
  BlockMap
  BlockMultiVector
  RowMatrixTransposer
  MatrixMatrix
  TypeStack
  inout
  PerformanceCGSolve
  )

ADD_SUBDIRECTORY(BugTests)

IF (Tpetra_ENABLE_Epetra)
  ADD_SUBDIRECTORY(EpetraRowMatrix)
ENDIF()

INCLUDE(AddSubdirectories)

# for TpetraTestingUtilities.hpp
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

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
  VbrMatrix
  RowMatrixTransposer
  MatrixMatrix
  TypeStack
  inout
  )

ADD_SUBDIRECTORY(BugTests)

IF (Tpetra_ENABLE_Epetra)
  ADD_SUBDIRECTORY(EpetraRowMatrix)
ENDIF()

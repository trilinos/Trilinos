INCLUDE(AddSubdirectories)

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

TRIBITS_ADD_LIBRARY(
  tpetratestlib
  SOURCES UnitTestMissingExplicitInstantiations.cpp
  DEPLIBS tpetra
  TESTONLY
)


IF (Tpetra_ENABLE_Epetra)
  ADD_SUBDIRECTORY(EpetraRowMatrix)
ENDIF()

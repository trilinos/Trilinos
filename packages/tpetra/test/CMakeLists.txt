INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  Platform
  Map
  Distributor
  Directory
  ImportExport
  ImportExport2
  MultiVector
  BasicPerfTest
  Utils
  CrsGraph
  CrsMatrix
  BlockMap
  BlockMultiVector
  VbrMatrix
  RowMatrixTransposer
  )

PACKAGE_ADD_LIBRARY(
  tpetratestlib
  SOURCES UnitTestMissingExplicitInstantiations.cpp
  DEPLIBS tpetra teuchos
  TESTONLY
)


IF (Tpetra_ENABLE_Epetra)
  ADD_SUBDIRECTORY(EpetraRowMatrix)
ENDIF()

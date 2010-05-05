INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  Platform
  BlockMap
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
  )

PACKAGE_ADD_LIBRARY(
  tpetratestlib
  SOURCES UnitTestMissingExplicitInstantiations.cpp
  TESTONLY
)


IF (Tpetra_ENABLE_Epetra)
  ADD_SUBDIRECTORY(EpetraRowMatrix)
ENDIF()

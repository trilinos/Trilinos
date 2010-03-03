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
  )

PACKAGE_ADD_LIBRARY(
  tpetratestlib
  SOURCES UnitTestMissingExplicitInstantiations.cpp
  TESTONLY
)


IF (Tpetra_ENABLE_Epetra AND Tpetra_HAVE_MPI)
  ADD_SUBDIRECTORY(EpetraRowMatrix)
ENDIF()


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
  BlockCrsGraph
  BlockMap
  BlockMultiVector
  VbrMatrix
  RowMatrixTransposer
  MatrixMatrix
  inout
  )

#IF (Tpetra_ENABLE_BUGTESTS)
ADD_SUBDIRECTORY(BugTests)
#ENDIF()

PACKAGE_ADD_LIBRARY(
  tpetratestlib
  SOURCES UnitTestMissingExplicitInstantiations.cpp
  DEPLIBS tpetra teuchos
  TESTONLY
)


IF (Tpetra_ENABLE_Epetra)
  ADD_SUBDIRECTORY(EpetraRowMatrix)
ENDIF()

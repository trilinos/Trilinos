INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  CompareWithAztecOO
  Container
  CrsRiluk
  IC
  ILU
  IlukGraph
  OverlappingRowMatrix
  PreconditionerFactory
  Relaxation
  SetParameters
  TestAll
  BlockCheby
  #UseTranspose
  )
IF (${PACKAGE_NAME}_USING_HYPRE AND Trilinos_ENABLE_Galeri AND
Trilinos_ENABLE_Teuchos AND Trilinos_ENABLE_AztecOO)
  ADD_SUBDIRECTORY(Hypre)
ENDIF()

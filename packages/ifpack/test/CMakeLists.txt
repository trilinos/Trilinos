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
IF (${PACKAGE_NAME}_USING_HYPRE)
  ADD_SUBDIRECTORY(Hypre)
ENDIF()

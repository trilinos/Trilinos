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
IF (HAVE_HYPRE)
  ADD_SUBDIRECTORY(Hypre)
ENDIF()

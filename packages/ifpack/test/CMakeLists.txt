INCLUDE(AddSubdirectories)

IF (NOT Trilinos_NO_32BIT_GLOBAL_INDICES)

ADD_SUBDIRECTORIES(
  CompareWithAztecOO
  Container
  CrsRiluk
  IC
  ILU
  IHSS-SORa
  IlukGraph
  OverlappingRowMatrix
  Performance
  PreconditionerFactory
  Relaxation
  SetParameters
  TriDi
  TestAll
  BlockCheby
  #UseTranspose
  SupportGraph
  )
IF (${PACKAGE_NAME}_ENABLE_AztecOO)
  ADD_SUBDIRECTORY(AztecOO)
ENDIF()
IF (HAVE_HYPRE)
  ADD_SUBDIRECTORY(Hypre)
ENDIF()

ENDIF()

IF (NOT Trilinos_NO_64BIT_GLOBAL_INDICES)

ADD_SUBDIRECTORIES(
  CompareWithAztecOO_LL
  Container_LL
  CrsRiluk_LL
  IC_LL
  ILU_LL
  IHSS-SORa_LL
  IlukGraph_LL
  OverlappingRowMatrix_LL
  PreconditionerFactory_LL
  Relaxation_LL
  SetParameters_LL
  TestAll_LL
  BlockCheby_LL
  )
IF (${PACKAGE_NAME}_ENABLE_AztecOO)
  ADD_SUBDIRECTORY(AztecOO_LL)
ENDIF()
IF (HAVE_HYPRE)
  ADD_SUBDIRECTORY(Hypre_LL)
ENDIF()

ENDIF()


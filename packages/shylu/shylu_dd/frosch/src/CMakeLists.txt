#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)


#SET(CURRENT_PACKAGE ShyLU)
#CONFIGURE_FILE("${Trilinos_SOURCE_DIR}/packages/Trilinos_DLLExportMacro.h.in"
#  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_DLLExportMacro.h)


#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")


# Set the base Trilinos directory in order to find Trilinos_version.h

ASSERT_DEFINED(PROJECT_SOURCE_DIR)
ASSERT_DEFINED(PROJECT_BINARY_DIR)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

#
# Include packages/src directory
#
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/CoarseSpaces)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/InterfaceSets)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SchwarzOperators)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SchwarzPreconditioners)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Tools)


APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

## AH 09/13/2017:
## There is a problem with the dependecies of the packages... We have to discuss this
#INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/amesos2/src)
#INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/amesos2/src/KLU2/Include)
#INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/amesos2/src/KLU2/Source)
#INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/amesos2/src)
#INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/muelu/adapters/belos)

APPEND_SET(HEADERS
  CoarseSpaces/FROSch_CoarseSpace_decl.hpp
  CoarseSpaces/FROSch_CoarseSpace_def.hpp
  CoarseSpaces/FROSch_GDSWInterfacePartitionOfUnity_decl.hpp
  CoarseSpaces/FROSch_GDSWInterfacePartitionOfUnity_def.hpp
  CoarseSpaces/FROSch_InterfacePartitionOfUnity_decl.hpp
  CoarseSpaces/FROSch_InterfacePartitionOfUnity_def.hpp
  CoarseSpaces/FROSch_LocalPartitionOfUnityBasis_decl.hpp
  CoarseSpaces/FROSch_LocalPartitionOfUnityBasis_def.hpp
  InterfaceSets/FROSch_DDInterface_decl.hpp
  InterfaceSets/FROSch_DDInterface_def.hpp
  InterfaceSets/FROSch_EntitySet_decl.hpp
  InterfaceSets/FROSch_EntitySet_def.hpp
  InterfaceSets/FROSch_InterfaceEntity_decl.hpp
  InterfaceSets/FROSch_InterfaceEntity_def.hpp
  SchwarzOperators/FROSch_AlgebraicOverlappingOperator_decl.hpp
  SchwarzOperators/FROSch_AlgebraicOverlappingOperator_def.hpp
  SchwarzOperators/FROSch_CoarseOperator_decl.hpp
  SchwarzOperators/FROSch_CoarseOperator_def.hpp
  SchwarzOperators/FROSch_GDSWCoarseOperator_decl.hpp
  SchwarzOperators/FROSch_GDSWCoarseOperator_def.hpp
  SchwarzOperators/FROSch_HarmonicCoarseOperator_decl.hpp
  SchwarzOperators/FROSch_HarmonicCoarseOperator_def.hpp
  SchwarzOperators/FROSch_IPOUHarmonicCoarseOperator_decl.hpp
  SchwarzOperators/FROSch_IPOUHarmonicCoarseOperator_def.hpp
  SchwarzOperators/FROSch_OverlappingOperator_decl.hpp
  SchwarzOperators/FROSch_OverlappingOperator_def.hpp
  SchwarzOperators/FROSch_RGDSWCoarseOperator_decl.hpp
  SchwarzOperators/FROSch_RGDSWCoarseOperator_def.hpp
  SchwarzOperators/FROSch_SchwarzOperator_decl.hpp
  SchwarzOperators/FROSch_SchwarzOperator_def.hpp
  SchwarzOperators/FROSch_SumOperator_decl.hpp
  SchwarzOperators/FROSch_SumOperator_def.hpp
  SchwarzPreconditioners/FROSch_AlgebraicOverlappingPreconditioner_decl.hpp
  SchwarzPreconditioners/FROSch_AlgebraicOverlappingPreconditioner_def.hpp
  SchwarzPreconditioners/FROSch_GDSWPreconditioner_decl.hpp
  SchwarzPreconditioners/FROSch_GDSWPreconditioner_def.hpp
  SchwarzPreconditioners/FROSch_OneLevelPreconditioner_decl.hpp
  SchwarzPreconditioners/FROSch_OneLevelPreconditioner_def.hpp
  SchwarzPreconditioners/FROSch_RGDSWPreconditioner_decl.hpp
  SchwarzPreconditioners/FROSch_RGDSWPreconditioner_def.hpp
  SchwarzPreconditioners/FROSch_SchwarzPreconditioner_decl.hpp
  SchwarzPreconditioners/FROSch_SchwarzPreconditioner_def.hpp
  SchwarzPreconditioners/FROSch_TwoLevelPreconditioner_decl.hpp
  SchwarzPreconditioners/FROSch_TwoLevelPreconditioner_def.hpp
  Tools/FROSch_ExtractSubmatrices_decl.hpp
  Tools/FROSch_ExtractSubmatrices_def.hpp
  Tools/FROSch_SubdomainSolver_decl.hpp
  Tools/FROSch_SubdomainSolver_def.hpp
  Tools/FROSch_Tools_decl.hpp
  Tools/FROSch_Tools_def.hpp
)

# Set sources

APPEND_SET(SOURCES
  dummy.cpp
)

TRIBITS_ADD_LIBRARY(
  shylu_ddfrosch
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )


INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/dofmngr)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/evaluators)

APPEND_SET(HEADERS
  Panzer_ArrayTraits.hpp
  Panzer_AssemblyEngine.hpp
  Panzer_AssemblyEngineT.hpp
  Panzer_AssemblyEngine_InArgs.hpp
  Panzer_AssemblyEngine_TemplateBuilder.hpp
  Panzer_AssemblyEngine_TemplateManager.hpp
  Panzer_Basis.hpp
  Panzer_BasisValues.hpp
  Panzer_BasisValuesT.hpp
  Panzer_BC.hpp
  Panzer_BCStrategy.hpp
  Panzer_BCStrategyT.hpp
  Panzer_BCStrategy_Base.hpp
  Panzer_BCStrategy_Factory.hpp
  Panzer_BCStrategy_Factory_Defines.hpp
  Panzer_BCStrategy_TemplateManager.hpp
  Panzer_CellData.hpp
  Panzer_ConfigDefs.hpp
  Panzer_Dimension.hpp
  Panzer_EquationSet.hpp
  Panzer_EquationSet_Base.hpp
  Panzer_EquationSet_DefaultImpl.hpp
  Panzer_EquationSet_DefaultImplT.hpp
  Panzer_EquationSet_Factory.hpp
  Panzer_EquationSet_Factory_Defines.hpp
  Panzer_EquationSet_TemplateManager.hpp
  Panzer_ExplicitTemplateInstantiation.hpp
  Panzer_FieldManagerBuilder.hpp
  Panzer_FieldManagerBuilderT.hpp
  Panzer_InputPhysicsBlock.hpp
  Panzer_IntegrationRule.hpp
  Panzer_IntegrationValues.hpp
  Panzer_IntegrationValuesT.hpp
  Panzer_IntrepidBasisFactory.hpp
  Panzer_MaterialModelEntry.hpp
  Panzer_ModelFactory.hpp
  Panzer_ModelFactory_TemplateManager.hpp
  Panzer_PhysicsBlock.hpp
  Panzer_Version.hpp
  Panzer_Workset.hpp
  Panzer_Workset_Builder.hpp
  Panzer_Workset_BuilderT.hpp
  dofmngr/Panzer_FieldPattern.hpp
  dofmngr/Panzer_IntrepidFieldPattern.hpp
  dofmngr/Panzer_GeometricAggFieldPattern.hpp
  dofmngr/Panzer_FieldAggPattern.hpp
  dofmngr/Panzer_ConnManager.hpp
  dofmngr/Panzer_DOFManager.hpp
  dofmngr/Panzer_DOFManagerT.hpp
  dofmngr/Panzer_UniqueGlobalIndexer.hpp
  dofmngr/Panzer_EpetraLinearObjFactory.hpp
  dofmngr/Panzer_EpetraLinearObjFactoryT.hpp
  evaluators/Panzer_Constant.hpp
  evaluators/Panzer_ConstantT.hpp
  evaluators/Panzer_GatherSolution_Epetra.hpp
  evaluators/Panzer_GatherSolution_EpetraT.hpp
  evaluators/Panzer_ScatterResidual_Epetra.hpp
  evaluators/Panzer_ScatterResidual_EpetraT.hpp
  evaluators/Panzer_ScatterDirichletResidual_Epetra.hpp
  evaluators/Panzer_ScatterDirichletResidual_EpetraT.hpp
  ) 

APPEND_SET(SOURCES
  Panzer_Basis.cpp
  Panzer_BC.cpp
  Panzer_BCStrategy.cpp
  Panzer_Dimension.cpp
  Panzer_EquationSet_DefaultImpl.cpp
  Panzer_InputPhysicsBlock.cpp
  Panzer_IntegrationRule.cpp
  Panzer_MaterialModelEntry.cpp
  Panzer_PhysicsBlock.cpp
  Panzer_Traits.cpp
  Panzer_Version.cpp
  dofmngr/Panzer_FieldPattern.cpp
  dofmngr/Panzer_IntrepidFieldPattern.cpp
  dofmngr/Panzer_GeometricAggFieldPattern.cpp
  dofmngr/Panzer_FieldAggPattern.cpp
  dofmngr/Panzer_DOFManager.cpp
  dofmngr/Panzer_EpetraLinearObjFactory.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  panzer
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )



INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

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
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/responses)

# auto generate headers, add to dependencies and include directories
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR})
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/dofmngr)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/evaluators)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/responses)
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_BINARY_DIR}/*.hpp)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/evaluators/*.hpp)
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/responses/*.hpp)

APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/evaluators/*.cpp)
APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/responses/*.cpp)

APPEND_SET(HEADERS
  Panzer_ArrayTraits.hpp
  Panzer_AssemblyEngine.hpp
  Panzer_AssemblyEngineT.hpp
  Panzer_AssemblyEngine_InArgs.hpp
  Panzer_AssemblyEngine_TemplateBuilder.hpp
  Panzer_AssemblyEngine_TemplateManager.hpp
  Panzer_Base.hpp
  Panzer_BasisIRLayout.hpp
  Panzer_PureBasis.hpp
  Panzer_BasisValues.hpp
  Panzer_BasisValuesT.hpp
  Panzer_BC.hpp
  Panzer_BCStrategy.hpp
  Panzer_BCStrategyT.hpp
  Panzer_BCStrategy_Base.hpp
  Panzer_BCStrategy_Factory.hpp
  Panzer_BCStrategy_Factory_Defines.hpp
  Panzer_BCStrategy_TemplateManager.hpp
  Panzer_BCStrategy_Dirichlet_DefaultImpl.hpp
  Panzer_BCStrategy_Dirichlet_DefaultImplT.hpp
  Panzer_BCStrategy_Dirichlet_Utilities.hpp
  Panzer_CellData.hpp
  Panzer_CloneableEvaluator.hpp
  Panzer_ClosureModel_Factory.hpp
  Panzer_ClosureModel_Factory_TemplateManager.hpp
  Panzer_ConfigDefs.hpp
  Panzer_Dimension.hpp
  Panzer_GlobalData.hpp
  Panzer_GlobalDataAcceptor.hpp
  Panzer_GlobalDataAcceptor_DefaultImpl.hpp
  Panzer_EquationSet.hpp
  Panzer_EquationSet_Base.hpp
  Panzer_EquationSet_DefaultImpl.hpp
  Panzer_EquationSet_DefaultImplT.hpp
  Panzer_EquationSet_Factory.hpp
  Panzer_EquationSet_Factory_Defines.hpp
  Panzer_EquationSet_TemplateManager.hpp
  Panzer_EvaluationTraits.hpp
  Panzer_ExplicitTemplateInstantiation.hpp
  Panzer_FieldManagerBuilder.hpp
  Panzer_FieldManagerBuilderT.hpp
  Panzer_FieldLibrary.hpp
  Panzer_InitialCondition_Builder.hpp
  Panzer_InitialCondition_BuilderT.hpp
  Panzer_InputEquationSet.hpp
  Panzer_InputPhysicsBlock.hpp
  Panzer_IntegrationRule.hpp
  Panzer_IntegrationValues.hpp
  Panzer_IntegrationValuesT.hpp
  Panzer_IntrepidBasisFactory.hpp
  Panzer_MaterialModelEntry.hpp
  Panzer_ModelEvaluator.hpp
  Panzer_ModelEvaluatorT.hpp
  Panzer_ModelEvaluator_Epetra.hpp
  Panzer_ModelEvaluator_Factory.hpp
  Panzer_ModelEvaluator_FactoryT.hpp
  Panzer_OutputStream.hpp
  Panzer_OutputStream_DefaultImpl.hpp
  Panzer_ParameterLibrary.hpp
  Panzer_ParameterLibraryAcceptor.hpp
  Panzer_ParameterLibraryAcceptor_DefaultImpl.hpp
  Panzer_ParameterLibraryUtilities.hpp
  Panzer_ParameterLibraryUtilitiesT.hpp
  Panzer_ParameterList_ObjectBuilders.hpp
  Panzer_PauseToAttach.hpp
  Panzer_PhysicsBlock.hpp
  Panzer_ScalarParameterEntry.hpp
  Panzer_Shards_Utilities.hpp
  Panzer_StlMap_Utilities.hpp
  Panzer_String_Utilities.hpp
  Panzer_Traits.hpp
  Panzer_Version.hpp
  Panzer_Workset.hpp
  Panzer_Workset_Builder.hpp
  Panzer_Workset_BuilderT.hpp
  Panzer_WorksetContainer.hpp
  Panzer_WorksetFactoryBase.hpp
  Panzer_Workset_Utilities.hpp
  Panzer_ZeroSensitivities.hpp
  dofmngr/Panzer_FieldPattern.hpp
  dofmngr/Panzer_IntrepidFieldPattern.hpp
  dofmngr/Panzer_GeometricAggFieldPattern.hpp
  dofmngr/Panzer_FieldAggPattern.hpp
  dofmngr/Panzer_ConnManager.hpp
  dofmngr/Panzer_DOFManager.hpp
  dofmngr/Panzer_DOFManagerT.hpp
  dofmngr/Panzer_BlockedDOFManager.hpp
  dofmngr/Panzer_BlockedDOFManagerT.hpp
  dofmngr/Panzer_UniqueGlobalIndexer_Utilities.hpp
  dofmngr/Panzer_UniqueGlobalIndexer_UtilitiesT.hpp
  dofmngr/Panzer_UniqueGlobalIndexer.hpp
  dofmngr/Panzer_EpetraLinearObjFactory.hpp
  dofmngr/Panzer_EpetraLinearObjFactoryT.hpp
  dofmngr/Panzer_EpetraLinearObjContainer.hpp
  dofmngr/Panzer_UniqueGlobalIndexerFactory.hpp
  dofmngr/Panzer_DOFManagerFactory.hpp
  dofmngr/Panzer_DOFManagerFactoryT.hpp
  dofmngr/Panzer_LinearObjFactory.hpp
  ) 

APPEND_SET(SOURCES
  Panzer_BasisIRLayout.cpp
  Panzer_PureBasis.cpp
  Panzer_BC.cpp
  Panzer_BCStrategy.cpp
  Panzer_Dimension.cpp
  Panzer_EquationSet_DefaultImpl.cpp
  Panzer_GlobalData.cpp
  Panzer_GlobalDataAcceptor_DefaultImpl.cpp
  Panzer_FieldLibrary.cpp
  Panzer_InputEquationSet.cpp
  Panzer_InputPhysicsBlock.cpp
  Panzer_IntegrationRule.cpp
  Panzer_MaterialModelEntry.cpp
  Panzer_ModelEvaluator_Epetra.cpp
  Panzer_OutputStream_DefaultImpl.cpp
  Panzer_ParameterList_ObjectBuilders.cpp
  Panzer_PhysicsBlock.cpp
  Panzer_String_Utilities.cpp
  Panzer_Traits.cpp
  Panzer_Version.cpp
  Panzer_WorksetContainer.cpp
  Panzer_Workset_Utilities.cpp
  dofmngr/Panzer_FieldPattern.cpp
  dofmngr/Panzer_IntrepidFieldPattern.cpp
  dofmngr/Panzer_GeometricAggFieldPattern.cpp
  dofmngr/Panzer_FieldAggPattern.cpp
  dofmngr/Panzer_DOFManager.cpp
  dofmngr/Panzer_BlockedDOFManager.cpp
  dofmngr/Panzer_EpetraLinearObjFactory.cpp
  )

IF(${PACKAGE_NAME}_ENABLE_Stokhos)
  APPEND_SET(HEADERS
    Panzer_SGUtilities.hpp
    dofmngr/Panzer_SGEpetraLinearObjFactory.hpp
    dofmngr/Panzer_SGEpetraLinearObjFactoryT.hpp
    dofmngr/Panzer_SGEpetraLinearObjContainer.hpp
    )
  APPEND_SET(SOURCES
    dofmngr/Panzer_SGEpetraLinearObjFactory.cpp
    dofmngr/Panzer_SGEpetraLinearObjContainer.cpp
    )
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  panzer
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )


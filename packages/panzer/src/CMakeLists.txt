
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
  Panzer_BCStrategy_Dirichlet_DefaultImpl.hpp
  Panzer_BCStrategy_Dirichlet_DefaultImplT.hpp
  Panzer_CellData.hpp
  Panzer_CloneableEvaluator.hpp
  Panzer_ClosureModel_Factory.hpp
  Panzer_ClosureModel_Factory_TemplateManager.hpp
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
  Panzer_ParameterList_ObjectBuilders.hpp
  Panzer_PauseToAttach.hpp
  Panzer_PhysicsBlock.hpp
  Panzer_ResponseLibrary.hpp
  Panzer_ResponseContainer.hpp
  Panzer_String_Utilities.hpp
  Panzer_Version.hpp
  Panzer_Workset.hpp
  Panzer_Workset_Builder.hpp
  Panzer_Workset_BuilderT.hpp
  Panzer_VectorTemplateManager.hpp
  Panzer_VectorTemplateIterator.hpp
  Panzer_ZeroSensitivities.hpp
  Panzer_AuxiliaryEvaluator_TemplateManager.hpp
  Panzer_AuxiliaryEvaluator_Factory.hpp
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
  dofmngr/Panzer_SGEpetraLinearObjFactory.hpp
  dofmngr/Panzer_SGEpetraLinearObjFactoryT.hpp
  dofmngr/Panzer_EpetraLinearObjContainer.hpp
  dofmngr/Panzer_SGEpetraLinearObjContainer.hpp
  dofmngr/Panzer_UniqueGlobalIndexerFactory.hpp
  dofmngr/Panzer_DOFManagerFactory.hpp
  dofmngr/Panzer_DOFManagerFactoryT.hpp
  dofmngr/Panzer_LinearObjFactory.hpp
  evaluators/Panzer_Constant.hpp
  evaluators/Panzer_ConstantT.hpp
  evaluators/Panzer_Dirichlet_Residual.hpp
  evaluators/Panzer_Dirichlet_ResidualT.hpp
  evaluators/Panzer_DOF.hpp
  evaluators/Panzer_DOFT.hpp
  evaluators/Panzer_DOFGradient.hpp
  evaluators/Panzer_DOFGradientT.hpp
  evaluators/Panzer_GatherSolution_Epetra.hpp
  evaluators/Panzer_GatherSolution_EpetraT.hpp
  evaluators/Panzer_GatherSolution_EpetraSG.hpp
  evaluators/Panzer_GatherSolution_EpetraSGT.hpp
  evaluators/Panzer_GlobalStatistics.hpp
  evaluators/Panzer_GlobalStatisticsT.hpp
  evaluators/Panzer_Integrator_BasisTimesScalar.hpp
  evaluators/Panzer_Integrator_BasisTimesScalarT.hpp
  evaluators/Panzer_Integrator_GradBasisDotVector.hpp
  evaluators/Panzer_Integrator_GradBasisDotVectorT.hpp
  evaluators/Panzer_Normals.hpp
  evaluators/Panzer_NormalsT.hpp
  evaluators/Panzer_ScalarToVector.hpp
  evaluators/Panzer_ScalarToVectorT.hpp
  evaluators/Panzer_ScatterResidual_Epetra.hpp
  evaluators/Panzer_ScatterResidual_EpetraT.hpp
  evaluators/Panzer_ScatterResidual_EpetraSG.hpp
  evaluators/Panzer_ScatterResidual_EpetraSGT.hpp
  evaluators/Panzer_ScatterResponses.hpp
  evaluators/Panzer_ScatterResponsesT.hpp
  evaluators/Panzer_ScatterDirichletResidual_Epetra.hpp
  evaluators/Panzer_ScatterDirichletResidual_EpetraT.hpp
  evaluators/Panzer_ScatterDirichletResidual_EpetraSG.hpp
  evaluators/Panzer_ScatterDirichletResidual_EpetraSGT.hpp
  evaluators/Panzer_ScatterInitialCondition_Epetra.hpp
  evaluators/Panzer_ScatterInitialCondition_EpetraT.hpp
  evaluators/Panzer_ScatterInitialCondition_EpetraSG.hpp
  evaluators/Panzer_ScatterInitialCondition_EpetraSGT.hpp
  evaluators/Panzer_TestScatter.hpp
  evaluators/Panzer_TestScatterT.hpp
  evaluators/Panzer_Sum.hpp
  evaluators/Panzer_SumT.hpp
  evaluators/Panzer_VectorToScalar.hpp
  evaluators/Panzer_VectorToScalarT.hpp
  ) 

APPEND_SET(SOURCES
  Panzer_Basis.cpp
  Panzer_BC.cpp
  Panzer_BCStrategy.cpp
  Panzer_Dimension.cpp
  Panzer_EquationSet_DefaultImpl.cpp
  Panzer_InputEquationSet.cpp
  Panzer_InputPhysicsBlock.cpp
  Panzer_IntegrationRule.cpp
  Panzer_MaterialModelEntry.cpp
  Panzer_ModelEvaluator_Epetra.cpp
  Panzer_ParameterList_ObjectBuilders.cpp
  Panzer_PhysicsBlock.cpp
  Panzer_ResponseLibrary.cpp
  Panzer_ResponseContainer.cpp
  Panzer_String_Utilities.cpp
  Panzer_Traits.cpp
  Panzer_Version.cpp
  dofmngr/Panzer_FieldPattern.cpp
  dofmngr/Panzer_IntrepidFieldPattern.cpp
  dofmngr/Panzer_GeometricAggFieldPattern.cpp
  dofmngr/Panzer_FieldAggPattern.cpp
  dofmngr/Panzer_DOFManager.cpp
  dofmngr/Panzer_BlockedDOFManager.cpp
  dofmngr/Panzer_EpetraLinearObjFactory.cpp
  dofmngr/Panzer_SGEpetraLinearObjFactory.cpp
  dofmngr/Panzer_SGEpetraLinearObjContainer.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  panzer
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )


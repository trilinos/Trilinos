# Only enable this example if PHX_ENABLE_DEVICE_DAG is enabled

IF(PHX_ENABLE_DEVICE_DAG)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING 
                            ${CMAKE_CURRENT_SOURCE_DIR})
  
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING 
                    ${CMAKE_CURRENT_SOURCE_DIR}/evaluators)

# Reuse files from other example
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING 
                            ${${PACKAGE_NAME}_SOURCE_DIR}/example/FiniteElementAssembly)

SET(example_fem_SOURCES
    ${${PACKAGE_NAME}_SOURCE_DIR}/example/FiniteElementAssembly/MyTraits.hpp
    ${${PACKAGE_NAME}_SOURCE_DIR}/example/FiniteElementAssembly/ExplicitTemplateInstantiation.hpp
    ${${PACKAGE_NAME}_SOURCE_DIR}/example/FiniteElementAssembly/Workset.hpp
    ${${PACKAGE_NAME}_SOURCE_DIR}/example/FiniteElementAssembly/WorksetBuilder.hpp
    ${${PACKAGE_NAME}_SOURCE_DIR}/example/FiniteElementAssembly/Mesh.hpp
    ${${PACKAGE_NAME}_SOURCE_DIR}/example/FiniteElementAssembly/Mesh.cpp
    ${${PACKAGE_NAME}_SOURCE_DIR}/example/FiniteElementAssembly/Dimension.hpp
    ${${PACKAGE_NAME}_SOURCE_DIR}/example/FiniteElementAssembly/LinearObjectFactory.hpp
    ${${PACKAGE_NAME}_SOURCE_DIR}/example/FiniteElementAssembly/CommandLineParser.hpp
    ${${PACKAGE_NAME}_SOURCE_DIR}/example/FiniteElementAssembly/PrintValues.hpp
    evaluators/Constant.hpp
    evaluators/Constant_Def.hpp
    evaluators/Constant.cpp
    evaluators/IntegrateDiffusionTerm.hpp
    evaluators/IntegrateDiffusionTerm_Def.hpp
    evaluators/IntegrateDiffusionTerm.cpp
    evaluators/IntegrateSourceTerm.hpp
    evaluators/IntegrateSourceTerm_Def.hpp
    evaluators/IntegrateSourceTerm.cpp
    evaluators/GatherSolution.hpp
    evaluators/GatherSolution_Def.hpp
    evaluators/GatherSolution.cpp
    evaluators/ProjectValueToQP.hpp
    evaluators/ProjectValueToQP_Def.hpp
    evaluators/ProjectValueToQP.cpp
    evaluators/ProjectGradientToQP.hpp
    evaluators/ProjectGradientToQP_Def.hpp
    evaluators/ProjectGradientToQP.cpp
    evaluators/ZeroContributedField.hpp
    evaluators/ZeroContributedField_Def.hpp
    evaluators/ZeroContributedField.cpp
    evaluators/ScatterResidual.hpp
    evaluators/ScatterResidual_Def.hpp
    evaluators/ScatterResidual.cpp
    )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_finite_element_assembly_device_dag
  SOURCES FiniteElementAssembly_DeviceDag.cpp ${example_fem_SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

ENDIF(PHX_ENABLE_DEVICE_DAG)

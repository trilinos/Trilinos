TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING 
                            ${CMAKE_CURRENT_SOURCE_DIR})
  
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING 
                    ${CMAKE_CURRENT_SOURCE_DIR}/evaluators)

SET(example_fem_SOURCES
    MyTraits.hpp
    ExplicitTemplateInstantiation.hpp
    Dimension.hpp
    Workset.hpp
    WorksetBuilder.hpp
    Mesh.hpp
    Mesh.cpp
    LinearObjectFactory.hpp
    CommandLineParser.hpp
    PrintValues.hpp
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

# new implementation
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_finite_element_assembly
  SOURCES FiniteElementAssembly.cpp ${example_fem_SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  example_fem_nonlinear_Jv_timings
#  SOURCES Example_FEM_Nonlinear_Jv_Timings.cpp 
#          ${example_fem_nonlinear_SOURCES}
#  COMM serial mpi
#  PASS_REGULAR_EXPRESSION "Run has completed successfully!"
#  )

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING 
                            ${CMAKE_CURRENT_SOURCE_DIR}/auxiliary_tests)
  
SET(example_element_linear_2D_SOURCES
    Mesh.hpp
    Mesh.cpp
    Workset.hpp
    WorksetBuilder.hpp
    )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_mesh_aux_test
  SOURCES auxiliary_tests/Example_FEM_Mesh_AuxTest.cpp ${example_element_linear_2D_SOURCES} ${PHALANX_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_workset_aux_test
  SOURCES     auxiliary_tests/Example_FEM_Workset_AuxTest.cpp ${example_element_linear_2D_SOURCES} ${PHALANX_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_lof_aux_test
  SOURCES     auxiliary_tests/Example_FEM_LinearObjectFactory_AuxTest.cpp ${example_element_linear_2D_SOURCES} ${PHALANX_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

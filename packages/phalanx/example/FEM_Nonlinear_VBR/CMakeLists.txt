
INCLUDE(TribitsAddExecutableAndTest)

ASSERT_DEFINED(Phalanx_ENABLE_Belos)
ASSERT_DEFINED(Phalanx_ENABLE_Epetra)
ASSERT_DEFINED(Phalanx_ENABLE_Ifpack)

IF (Phalanx_ENABLE_Epetra AND Phalanx_ENABLE_Belos AND Phalanx_ENABLE_Ifpack AND Phalanx_ENABLE_ML AND Phalanx_ENABLE_Amesos)

  TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
  
  TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/evaluators)
  
  TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../FEM_Nonlinear)
  
  SET(example_fem_nonlinear_vbr_SOURCES
    Example_FEM_Nonlinear_VBR.cpp
    LinearObjectFactoryVBR.hpp
    LinearObjectFactoryVBR.cpp
    ../FEM_Nonlinear/Traits.hpp
    ../FEM_Nonlinear/Traits.cpp
    ../FEM_Nonlinear/FactoryTraits.hpp
    ../FEM_Nonlinear/Element_Linear2D.hpp
    ../FEM_Nonlinear/Element_Linear2D.cpp
    ../FEM_Nonlinear/ExplicitTemplateInstantiation.hpp
    ../FEM_Nonlinear/Dimension.hpp
    ../FEM_Nonlinear/Dimension.cpp
    ../FEM_Nonlinear/MeshBuilder.hpp
    ../FEM_Nonlinear/MeshBuilder.cpp
    ../FEM_Nonlinear/Workset.hpp
    evaluators/Evaluator_Constant.hpp
    evaluators/Evaluator_Constant_Def.hpp
    evaluators/Evaluator_Constant.cpp
    evaluators/Evaluator_Equations.hpp
    evaluators/Evaluator_Equations_Def.hpp
    evaluators/Evaluator_Equations.cpp
    evaluators/Evaluator_FEInterpolation.hpp
    evaluators/Evaluator_FEInterpolation_Def.hpp
    evaluators/Evaluator_FEInterpolation.cpp
    evaluators/Evaluator_GatherSolution.hpp
    evaluators/Evaluator_GatherSolution_Def.hpp
    evaluators/Evaluator_GatherSolution.cpp
    evaluators/Evaluator_ScatterResidual.hpp
    evaluators/Evaluator_ScatterResidual_Def.hpp
    evaluators/Evaluator_ScatterResidual.cpp
    )
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    example_fem_nonlinear_vbr
    SOURCES ${example_fem_nonlinear_vbr_SOURCES}
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Run has completed successfully!"
    )

  TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/auxiliary_tests)
  
  SET(example_epetra_vbr_test_SOURCES
    auxiliary_tests/Example_Epetra_VBR_Test.cpp
    ../FEM_Nonlinear/Traits.hpp
    ../FEM_Nonlinear/Traits.cpp
    ../FEM_Nonlinear/FactoryTraits.hpp
    ../FEM_Nonlinear/Element_Linear2D.hpp
    ../FEM_Nonlinear/Element_Linear2D.cpp
    ../FEM_Nonlinear/ExplicitTemplateInstantiation.hpp
    ../FEM_Nonlinear/Dimension.hpp
    ../FEM_Nonlinear/Dimension.cpp
    ../FEM_Nonlinear/Workset.hpp
    evaluators/Evaluator_Constant.hpp
    evaluators/Evaluator_Constant_Def.hpp
    evaluators/Evaluator_Constant.cpp
    evaluators/Evaluator_Equations.hpp
    evaluators/Evaluator_Equations_Def.hpp
    evaluators/Evaluator_Equations.cpp
    evaluators/Evaluator_FEInterpolation.hpp
    evaluators/Evaluator_FEInterpolation_Def.hpp
    evaluators/Evaluator_FEInterpolation.cpp
    evaluators/Evaluator_GatherSolution.hpp
    evaluators/Evaluator_GatherSolution_Def.hpp
    evaluators/Evaluator_GatherSolution.cpp
    evaluators/Evaluator_ScatterResidual.hpp
    evaluators/Evaluator_ScatterResidual_Def.hpp
    evaluators/Evaluator_ScatterResidual.cpp
    )
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    example_epetra_vbr_test
    SOURCES ${example_epetra_vbr_test_SOURCES}
    COMM serial mpi
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Run has completed successfully!"
    )
  
ENDIF(Phalanx_ENABLE_Epetra AND Phalanx_ENABLE_Belos AND Phalanx_ENABLE_Ifpack AND Phalanx_ENABLE_ML AND Phalanx_ENABLE_Amesos)

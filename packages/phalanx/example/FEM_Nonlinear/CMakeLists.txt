
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(Phalanx_ENABLE_Belos)
ASSERT_DEFINED(Phalanx_ENABLE_Epetra)
ASSERT_DEFINED(Phalanx_ENABLE_Ifpack)

IF (Phalanx_ENABLE_Epetra AND Phalanx_ENABLE_Belos AND Phalanx_ENABLE_Ifpack AND Phalanx_ENABLE_ML)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/evaluators)
  
  SET(example_fem_nonlinear_SOURCES
    Example_FEM_Nonlinear.cpp
    Traits.hpp
    Traits.cpp
    FactoryTraits.hpp
    Element_Linear2D.hpp
    Element_Linear2D.cpp
    ExplicitTemplateInstantiation.hpp
    Dimension.hpp
    Dimension.cpp
    Workset.hpp
    MeshBuilder.hpp
    MeshBuilder.cpp
    LinearObjectFactory.hpp
    LinearObjectFactory.cpp
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
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    example_fem_nonlinear
    SOURCES ${example_fem_nonlinear_SOURCES}
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Run has completed successfully!"
    )

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/auxiliary_tests)
  
  SET(example_element_linear_2D_SOURCES
    auxiliary_tests/Example_Element_Linear_2D.cpp
    Dimension.hpp
    Dimension.cpp
    Element_Linear2D.hpp
    Element_Linear2D.cpp
    MeshBuilder.hpp
    MeshBuilder.cpp
    )
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    example_element_linear_2D
    SOURCES ${example_element_linear_2D_SOURCES}
    COMM serial mpi
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Run has completed successfully!"
    )
  
ENDIF(Phalanx_ENABLE_Epetra AND Phalanx_ENABLE_Belos AND Phalanx_ENABLE_Ifpack AND Phalanx_ENABLE_ML)

INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/evaluators)

SET(example_SOURCES
  Example.cpp
  AlgebraicTypes.hpp
  Traits.hpp
  Traits.cpp
  FactoryTraits.hpp
  Cell.hpp
  Cell.cpp
  Workset.hpp
  evaluators/Evaluator_Constant.hpp
  evaluators/Evaluator_Constant_Def.hpp
  evaluators/Evaluator_Density.hpp
  evaluators/Evaluator_Density_Def.hpp
  evaluators/Evaluator_EnergyFlux_Fourier.hpp
  evaluators/Evaluator_EnergyFlux_Fourier_Def.hpp
  evaluators/Evaluator_FEInterpolation.hpp
  evaluators/Evaluator_FEInterpolation_Def.hpp
  evaluators/Evaluator_NonlinearSource.hpp
  evaluators/Evaluator_NonlinearSource_Def.hpp
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  example
  SOURCES ${example_SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Run has completed successfully!"
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/auxiliary_tests)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  example_debug_names
  SOURCES auxiliary_tests/Example_DebugNames.cpp 
          AlgebraicTypes.hpp 
          Traits.hpp 
          Traits.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Run has completed successfully!"
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  example_algebraic_types
  SOURCES auxiliary_tests/Example_AlgebraicTypes.cpp
          AlgebraicTypes.hpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Run has completed successfully!"
  )

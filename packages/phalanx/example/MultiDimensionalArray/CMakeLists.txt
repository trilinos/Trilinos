
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/evaluators)

SET(example_SOURCES
  Example_MultiDimensionalArray.cpp
  Traits.hpp
  Traits.cpp
  FactoryTraits.hpp
  Cell.hpp
  Cell.cpp
  Dimension.hpp
  Dimension.cpp
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

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_multidimensional_array
  SOURCES ${example_SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Run has completed successfully!"
  )

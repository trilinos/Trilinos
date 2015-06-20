

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/evaluators)


SET(example_SOURCES
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

SET(example_DFAD_SOURCES Example_DFAD.cpp TraitsDFAD.hpp ${example_SOURCES})
SET(example_SFAD_SOURCES Example_SFAD.cpp TraitsSFAD.hpp ${example_SOURCES})
SET(example_SLFAD_SOURCES Example_SLFAD.cpp TraitsSLFAD.hpp ${example_SOURCES})


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_dfad
  SOURCES ${example_DFAD_SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Run has completed successfully!"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_sfad
  SOURCES ${example_SFAD_SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Run has completed successfully!"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_slfad
  SOURCES ${example_SLFAD_SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Run has completed successfully!"
  )

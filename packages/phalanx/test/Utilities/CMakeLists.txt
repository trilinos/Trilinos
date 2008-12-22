
INCLUDE(TrilinosPackageLibraryHelpers)
INCLUDE(Trilinos_Add_Executable_And_Test)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  CellData.hpp
  AlgebraicTypes.hpp
  Traits.hpp
  Traits_ContiguousAllocator.hpp
  FactoryTraits.hpp
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

APPEND_SET(SOURCES
  CellData.cpp
  Traits.cpp
  Traits_ContiguousAllocator.cpp
  )

TRILINOS_PACKAGE_ADD_LIBRARY(
  phalanx_test_utilities
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  TESTONLY
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  test_utilitites_algebraic_types
  SOURCES Test_Utilities_AlgebraicTypes.cpp
  DEPLIBS phalanx_test_utilities
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Test passed!"
  )
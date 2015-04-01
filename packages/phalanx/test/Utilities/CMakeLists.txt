INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddExecutableAndTest)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})


APPEND_SET(HEADERS
  CellData.hpp
  Traits.hpp
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
  )

#ENDIF ()

TRIBITS_ADD_LIBRARY(
  phalanx_test_utilities
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  TESTONLY
  )


INCLUDE(PackageLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  OptiPack_DiagonalQuadraticResponseOnlyModelEvaluator.hpp
  OptiPack_DiagonalQuadraticResponseOnlyModelEvaluator_decl.hpp
  OptiPack_DiagonalQuadraticResponseOnlyModelEvaluator_def.hpp
  )

SET(SOURCES ${SOURCES}
  OptiPack_DiagonalQuadraticResponseOnlyModelEvaluator.cpp
  )

PACKAGE_ADD_LIBRARY(
  optipack_test_models
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  TESTONLY
  )

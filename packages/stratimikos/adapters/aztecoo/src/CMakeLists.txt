
INCLUDE(TribitsLibraryMacros)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS
  AztecOOParameterList.hpp
  Thyra_AztecOOLinearOpWithSolve.hpp
  Thyra_AztecOOLinearOpWithSolveFactory.hpp
  )

SET(SOURCES
  AztecOOParameterList.cpp
  Thyra_AztecOOLinearOpWithSolve.cpp
  Thyra_AztecOOLinearOpWithSolveFactory.cpp
  )

TRIBITS_ADD_LIBRARY(
  stratimikosaztecoo
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

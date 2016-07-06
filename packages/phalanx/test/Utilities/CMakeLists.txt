INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddExecutableAndTest)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})


APPEND_SET(HEADERS
  CellData.hpp
  Traits.hpp
  Evaluator_MockDAG.hpp
  Evaluator_MockDAG_Def.hpp
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

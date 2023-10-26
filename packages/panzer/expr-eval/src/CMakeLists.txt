SET(HEADERS
    Panzer_ExprEval.hpp
    Panzer_ExprEval_impl.hpp
   )
SET(SOURCES
    Panzer_ExprEval.cpp
   )

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_LIBRARY(
  panzer-expr-eval
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

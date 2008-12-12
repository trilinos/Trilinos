
INCLUDE(TrilinosPackageLibraryHelpers)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS
  BelosThyraAdapter.hpp
  Thyra_BelosLinearOpWithSolveDecl.hpp
  Thyra_BelosLinearOpWithSolveFactoryDecl.hpp
  Thyra_BelosLinearOpWithSolveFactory.hpp
  Thyra_BelosLinearOpWithSolve.hpp 
  )

SET(SOURCES
  Thyra_BelosAdapters.cpp
  )

TRILINOS_PACKAGE_ADD_LIBRARY(
  stratimikosbelos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )


INCLUDE(PackageLibraryMacros)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS
  BelosThyraAdapter.hpp
  Thyra_BelosLinearOpWithSolve.hpp
  Thyra_BelosLinearOpWithSolve_decl.hpp
  Thyra_BelosLinearOpWithSolve_def.hpp
  Thyra_BelosLinearOpWithSolveFactory.hpp
  Thyra_BelosLinearOpWithSolveFactory_decl.hpp
  Thyra_BelosLinearOpWithSolveFactory_def.hpp
  )

SET(SOURCES
  Thyra_BelosAdapters.cpp
  Thyra_BelosLinearOpWithSolve.cpp
  Thyra_BelosLinearOpWithSolveFactory.cpp
  )

PACKAGE_ADD_LIBRARY(
  stratimikosbelos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

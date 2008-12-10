
INCLUDE(TrilinosPackageLibraryHelpers)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS
  Thyra_AmesosTypes.hpp
  Thyra_AmesosLinearOpWithSolve.hpp
  Thyra_AmesosLinearOpWithSolveFactory.hpp
  )

SET(SOURCES
  Thyra_AmesosTypes.cpp
  Thyra_AmesosLinearOpWithSolve.cpp
  Thyra_AmesosLinearOpWithSolveFactory.cpp
  )

TRILINOS_PACKAGE_ADD_LIBRARY(
  stratimikosamesos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

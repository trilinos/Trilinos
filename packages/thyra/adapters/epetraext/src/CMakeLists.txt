
INCLUDE(TrilinosPackageLibraryHelpers)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/model_evaluator)

SET(HEADERS
  model_evaluator/Thyra_EpetraModelEvaluator.hpp
  )

SET(SOURCES
  model_evaluator/Thyra_EpetraModelEvaluator.cpp
  )

TRILINOS_PACKAGE_ADD_LIBRARY(
  thyraepetraext
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS thyraepetra
  )

TRILINOS_PACKAGE_EXPORT_DEPENDENCY_VARIABLES()

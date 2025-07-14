TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src-epetra)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../src-epetra)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../headers)


APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
)

APPEND_SET(HEADERS
  ../headers/Galeri_ConfigDefs.h
  ../headers/Galeri_Exception.h

  Galeri_XpetraParameters.hpp

  Galeri_XpetraProblemFactory.hpp
  Galeri_XpetraProblemFactory_Helmholtz.hpp
  Galeri_XpetraMatrixTypes.hpp
  Galeri_XpetraMatrixTypes_Helmholtz.hpp
  Galeri_StencilProblems.hpp
  Galeri_StencilProblems_Helmholtz.hpp
  Galeri_Elasticity2DProblem.hpp
  Galeri_Elasticity3DProblem.hpp
  Galeri_Problem_Helmholtz.hpp
  Galeri_HelmholtzFEM2DProblem.hpp
  Galeri_HelmholtzFEM3DProblem.hpp
  Galeri_VelocityModel.hpp

  Galeri_VectorTraits.hpp
  Galeri_MatrixTraits.hpp
  Galeri_MapTraits.hpp
  Galeri_MultiVectorTraits.hpp

  Galeri_XpetraUtils.hpp
  Galeri_Problem.hpp

  Galeri_XpetraMaps.hpp
  Galeri_XpetraCartesian.hpp
  )

APPEND_SET(SOURCES
  mac.cpp
  )

TRIBITS_ADD_LIBRARY(
  galeri-xpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INCLUDE(TribitsLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/Utils)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/Headers)

APPEND_SET(HEADERS
  Galeri_XpetraParameters.hpp

  Galeri_XpetraProblemFactory.hpp
  Galeri_XpetraMatrixTypes.hpp
  Galeri_StencilProblems.hpp
  Galeri_Elasticity2DProblem.hpp
  Galeri_Elasticity3DProblem.hpp

  Galeri_VectorTraits.hpp
  Galeri_MatrixTraits.hpp
  Galeri_MapTraits.hpp
  Galeri_MultiVectorTraits.hpp
  
  Galeri_XpetraUtils.hpp
  Galeri_Problem.hpp

  Galeri_XpetraMaps.hpp
  Maps/Galeri_XpetraCartesian1D.hpp
  Maps/Galeri_XpetraCartesian2D.hpp
  Maps/Galeri_XpetraCartesian3D.hpp
  )

APPEND_SET(SOURCES
  mac.cpp
  )

TRIBITS_ADD_LIBRARY(
  galeri-xpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

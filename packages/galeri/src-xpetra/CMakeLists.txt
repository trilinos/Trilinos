INCLUDE(TribitsLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/Utils)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/Headers)

APPEND_SET(HEADERS
  Galeri_XpetraParameters.hpp
  Galeri_XpetraMatrixFactory.hpp  
  Galeri_XpetraMatrixTypes.hpp
  Galeri_XpetraVectorTraits.hpp
  )

APPEND_SET(SOURCES
  mac.cpp
  )

TRIBITS_ADD_LIBRARY(
  galeri-xpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INCLUDE(TribitsLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/Utils)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/Headers)

APPEND_SET(HEADERS
  MueLu_GalleryParameters.hpp
  MueLu_MatrixFactory.hpp  
  MueLu_MatrixTypes.hpp
  MueLu_VectorTraits.hpp
  )

APPEND_SET(SOURCES
  mac.cpp
  )

TRIBITS_ADD_LIBRARY(
  muelu-gallery
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS muelu # for timers etc.
  )

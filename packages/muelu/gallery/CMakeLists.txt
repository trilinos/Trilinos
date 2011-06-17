
INCLUDE(PackageLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../utils)

APPEND_SET(HEADERS
  MueLu_GalleryParameters.hpp
  MueLu_MatrixFactory.hpp  
  MueLu_MatrixTypes.hpp
  )

APPEND_SET(SOURCES
  mac.cpp
  )

PACKAGE_ADD_LIBRARY(
  muelu-gallery
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

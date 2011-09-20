INCLUDE(PackageLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS 
  belos/BelosMueLuAdapter.hpp 
  belos/BelosMueLuAdapterMultiVector.hpp
  )

#APPEND_SET(SOURCES mac.cpp)

PACKAGE_ADD_LIBRARY(
  muelu-adapters
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

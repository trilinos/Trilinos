INCLUDE(PackageLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/belos)

APPEND_SET(HEADERS
  belos/BelosOperatorT.hpp
  belos/BelosXpetraAdapter.hpp
  belos/BelosXpetraAdapterMultiVector.hpp
  belos/BelosXpetraAdapterOperator.hpp
  belos/BelosMueLuAdapter.hpp
  )

APPEND_SET(SOURCES mac.cpp)

PACKAGE_ADD_LIBRARY(
  muelu-adapters
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

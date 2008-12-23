
INCLUDE(PackageLibraryMacros)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS
  Thyra_MLPreconditionerFactory.hpp
  )

SET(SOURCES
  Thyra_MLPreconditionerFactory.cpp
  )

PACKAGE_ADD_LIBRARY(
  stratimikosml
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

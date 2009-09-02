
INCLUDE(PackageLibraryMacros)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS
  Thyra_IfpackPreconditionerFactory.hpp
  )

SET(SOURCES
  Thyra_IfpackPreconditionerFactory.cpp
  )

PACKAGE_ADD_LIBRARY(
  stratimikosifpack
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

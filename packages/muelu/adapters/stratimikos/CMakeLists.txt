INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS
  Thyra_MueLuPreconditionerFactory.hpp
  )

SET(SOURCES
  Thyra_MueLuPreconditionerFactory.cpp
  )

TRIBITS_ADD_LIBRARY(
  stratimikosmuelu
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS muelu
  )

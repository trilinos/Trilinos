set(HEADERS "")
set(SOURCES "")

tribits_configure_file(${PACKAGE_NAME}_config.h)
tribits_include_directories(${CMAKE_CURRENT_BINARY_DIR})

tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

append_set(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  Package3.hpp
  )
append_set(SOURCES
  Package3.cpp
  )

tribits_add_library(
  package3
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

tribits_add_executable(package3-prg  NOEXEPREFIX  NOEXESUFFIX
  SOURCES  Package3_Prg.cpp
  INSTALLABLE
  )

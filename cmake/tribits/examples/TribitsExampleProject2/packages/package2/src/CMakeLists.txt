set(HEADERS "")
set(SOURCES "")

tribits_configure_file(${PACKAGE_NAME}_config.h)
tribits_include_directories(${CMAKE_CURRENT_BINARY_DIR})

tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

append_set(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  Package2.hpp
  )
append_set(SOURCES
  Package2.cpp
  )

tribits_add_library(
  package2
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

tribits_add_executable(package2-prg  NOEXEPREFIX  NOEXESUFFIX
  SOURCES  Package2_Prg.cpp
  INSTALLABLE
  )

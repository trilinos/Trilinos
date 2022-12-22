set(HEADERS "")
set(SOURCES "")

tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

append_set(HEADERS
  Package1.hpp
  )
append_set(SOURCES
  Package1.cpp
  )

tribits_add_library(
  package1
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

tribits_add_executable(package1-prg  NOEXEPREFIX  NOEXESUFFIX
  SOURCES  Package1_Prg.cpp
  INSTALLABLE
  )

tribits_configure_file(${PACKAGE_NAME}_config.h)

add_subdirectory(vector)

set(HEADERS "")
set(SOURCES "")

tribits_include_directories(${CMAKE_CURRENT_BINARY_DIR})
append_set(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})
append_set(HEADERS
  MixedLang.hpp
  Ray_Tracer.hh
  Ray.hh
  )
append_set(SOURCES
  MixedLang.cpp
  Ray_Tracer.cc
  Parameters.f90
  Ray_Tracer_Kernel.f90
  )

tribits_add_library(
  mixedlang
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS mixedlang_vector
  )

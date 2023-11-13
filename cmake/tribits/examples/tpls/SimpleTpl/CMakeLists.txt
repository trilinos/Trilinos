cmake_minimum_required(VERSION 3.23.0)

project(SimpleTpl
  DESCRIPTION
    "Simple TPL with a header file and a library"
  LANGUAGES C CXX
  )

add_library(simpletpl SimpleTpl.hpp SimpleTpl.cpp)
target_include_directories(simpletpl
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )

include(GNUInstallDirs)

install(
  TARGETS simpletpl
  EXPORT SimpleTpl
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

install(
  FILES SimpleTpl.hpp
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )

install(EXPORT SimpleTpl
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SimpleTpl
  NAMESPACE SimpleTpl::
  FILE SimplTplConfig.cmake
  )

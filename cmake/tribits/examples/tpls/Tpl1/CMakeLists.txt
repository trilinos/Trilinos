cmake_minimum_required(VERSION 3.23.0)

project(Tpl1
  DESCRIPTION
    "Simple TPL with one header file, one library, and no deps"
  LANGUAGES C CXX
  )

add_library(tpl1 Tpl1.hpp Tpl1.cpp)
target_include_directories(tpl1
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )

include(GNUInstallDirs)

install(
  TARGETS tpl1
  EXPORT ${PROJECT_NAME}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

install(
  FILES Tpl1.hpp
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )

install(EXPORT ${PROJECT_NAME}
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  NAMESPACE tpl1::
  FILE ${PROJECT_NAME}Config.cmake
  )

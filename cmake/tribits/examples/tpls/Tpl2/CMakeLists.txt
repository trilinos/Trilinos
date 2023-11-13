cmake_minimum_required(VERSION 3.23.0)

project(Tpl2
  DESCRIPTION
    "Simple TPL with two header file, two libraries, and depends on Tpl1"
  LANGUAGES C CXX
  )

find_package(Tpl1 REQUIRED)


add_library(tpl2a Tpl2a.hpp Tpl2a.cpp)
target_include_directories(tpl2a
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(tpl2a PRIVATE tpl1::tpl1)

add_library(tpl2b Tpl2b.hpp Tpl2b.cpp)
target_include_directories(tpl2b
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# NOTE: The above libs tpl2a and tpl2b do not depend on each other.  They just
# depend on tpl1::tpl1.

include(GNUInstallDirs)

install(
  TARGETS tpl2a tpl2b
  EXPORT ${PROJECT_NAME}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

install(
  FILES Tpl2a.hpp Tpl2b.hpp
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )

install(EXPORT ${PROJECT_NAME}
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  NAMESPACE tpl2::
  FILE ${PROJECT_NAME}ConfigTargets.cmake
  )

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Tpl2Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Tpl2Config.install.cmake"
  @ONLY
  )

install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Tpl2Config.install.cmake"
  RENAME "Tpl2Config.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  )

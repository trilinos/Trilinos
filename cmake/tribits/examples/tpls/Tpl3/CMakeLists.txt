cmake_minimum_required(VERSION 3.23.0)

project(Tpl3
  DESCRIPTION
    "Simple TPL with one header file, one library, and depends on Tpl2"
  LANGUAGES C CXX
  )

find_package(Tpl2 REQUIRED)

add_library(tpl3 Tpl3.hpp Tpl3.cpp)
target_include_directories(tpl3
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(tpl3 PRIVATE tpl2::tpl2a tpl2::tpl2b)
# Does CMake maintain the above order?

include(GNUInstallDirs)

install(
  TARGETS tpl3
  EXPORT ${PROJECT_NAME}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

install(
  FILES Tpl3.hpp
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )

install(EXPORT ${PROJECT_NAME}
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  NAMESPACE tpl3::
  FILE ${PROJECT_NAME}ConfigTargets.cmake
  )

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Tpl3Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Tpl3Config.install.cmake"
  @ONLY
  )

install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Tpl3Config.install.cmake"
  RENAME "Tpl3Config.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  )

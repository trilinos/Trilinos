cmake_minimum_required(VERSION 3.23.0)

project(Tpl4
  DESCRIPTION
    "Simple TPL with one header file, no libraries, and depends on Tpl2 and Tpl3"
  LANGUAGES C CXX
  )

find_package(Tpl2 REQUIRED)
find_package(Tpl3 REQUIRED)

add_library(tpl4 INTERFACE)
target_include_directories(tpl4
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(tpl4
  INTERFACE tpl3::tpl3
  INTERFACE tpl2::tpl2a tpl2::tpl2b
  )
# Does CMake maintain the above order?

include(GNUInstallDirs)

install(
  TARGETS tpl4
  EXPORT ${PROJECT_NAME}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

install(
  FILES Tpl4.hpp
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )

install(EXPORT ${PROJECT_NAME}
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  NAMESPACE tpl4::
  FILE ${PROJECT_NAME}ConfigTargets.cmake
  )

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Tpl4Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Tpl4Config.install.cmake"
  @ONLY
  )

install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Tpl4Config.install.cmake"
  RENAME "Tpl4Config.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  )

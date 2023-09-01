cmake_minimum_required(VERSION 3.23.0 FATAL_ERROR)

macro(include_raw_cmake_build)
  if (NOT COMMAND tribits_project)
    include("${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.raw.cmake"
       NO_POLICY_SCOPE)
    return()
  endif()
endmacro()

if (COMMAND tribits_package)

  # Being processed as a TriBITS package
  tribits_package(Package1)
  add_subdirectory(src)
  tribits_add_test_directories(test)
  tribits_package_postprocess()

else()

  message("Configuring raw CMake project Package1")
  project(Package1 LANGUAGES C CXX)
  include(GNUInstallDirs)
  find_package(Tpl1 CONFIG REQUIRED)
  add_subdirectory(src)
  if (Package1_ENABLE_TESTS)
    include(CTest)
    add_subdirectory(test)
  endif()

  # Generate the all_libs target
  add_library(Package1_all_libs INTERFACE)
  set_target_properties(Package1_all_libs PROPERTIES
    EXPORT_NAME all_libs)
  target_link_libraries(Package1_all_libs INTERFACE Package1_package1)
  install(
    TARGETS Package1_all_libs
    EXPORT ${PROJECT_NAME}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )

  # Generate and install the Package1Config.cmake file
  install(EXPORT ${PROJECT_NAME}
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
    NAMESPACE ${PROJECT_NAME}::
    FILE ${PROJECT_NAME}ConfigTargets.cmake )
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/raw/Package1Config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Package1Config.install.cmake"
    @ONLY )
  install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Package1Config.install.cmake"
    RENAME "Package1Config.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME} )

endif()



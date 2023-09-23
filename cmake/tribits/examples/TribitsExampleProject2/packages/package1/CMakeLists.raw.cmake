cmake_minimum_required(VERSION 3.23.0 FATAL_ERROR)

if (COMMAND tribits_package)
  message("Configuring raw CMake package Package1")
else()
  message("Configuring raw CMake project Package1")
endif()

# Standard project-level stuff
project(Package1 LANGUAGES C CXX)
include(GNUInstallDirs)
find_package(Tpl1 CONFIG REQUIRED)
add_subdirectory(src)
if (Package1_ENABLE_TESTS)
  include(CTest)
  include("${CMAKE_CURRENT_LIST_DIR}/cmake/raw/EnableTribitsTestSupport.cmake")
  add_subdirectory(test)
endif()

# Stuff that TriBITS does automatically
include("${CMAKE_CURRENT_LIST_DIR}/cmake/raw/DefineAllLibsTarget.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cmake/raw/GeneratePackageConfigFileForBuildDir.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cmake/raw/GeneratePackageConfigFileForInstallDir.cmake")

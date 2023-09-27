cmake_minimum_required(VERSION 3.23.0 FATAL_ERROR)

set(Package1_USE_RAW_CMAKE  OFF  CACHE  BOOL
  "Use raw CMake for package build, even if TriBITS could be used.")

# Macro to select the TriBITS or the raw CMake build system
macro(include_cmakelists_file)
  if ((NOT COMMAND tribits_project) OR Package1_USE_RAW_CMAKE)
    include("${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.raw.cmake"
       NO_POLICY_SCOPE)
  else()
    include("${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.tribits.cmake"
       NO_POLICY_SCOPE)
  endif()
endmacro()

# Pull in the base CMakeLists.txt file variant
include_cmakelists_file()

cmake_minimum_required(VERSION 3.10.0 FATAL_ERROR)
include(${CMAKE_CURRENT_SOURCE_DIR}/ProjectName.cmake)

# Called at the top of every CMakeLists.txt file
macro(include_tribits_build)
  if (${PROJECT_NAME}_TRIBITS_DIR)
    include("${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.tribits.cmake")
    return()
  endif()
endmacro()

if (${PROJECT_NAME}_TRIBITS_DIR)

  # TriBITS CMake project
  project(${PROJECT_NAME} NONE)
  include("${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake")
  # Only one package in this simple project so just enable it :-)
  set(${PROJECT_NAME}_ENABLE_HelloWorld ON CACHE BOOL "" FORCE)
  tribits_project()

else()

  # Raw CMake project
  project(RawHelloWorld)
  enable_testing()
  add_subdirectory(hello_world)

endif()

# NOTE: The cmake_minimum_required() and project() commands must be executed
# in the base CMakeLists.txt file and *NOT* in an included 

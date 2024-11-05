# Example of a simple project that uses TriBITS
cmake_minimum_required(VERSION 3.23.0 FATAL_ERROR)
include(${CMAKE_CURRENT_SOURCE_DIR}/ProjectName.cmake)
project(${PROJECT_NAME} NONE)
set(${PROJECT_NAME}_TRIBITS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../..
  CACHE PATH "TriBITS base directory (default assumes in TriBITS source tree).")
include("${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake")
# Only one package in this simple project so just enable it :-)
set(${PROJECT_NAME}_ENABLE_HelloWorld ON CACHE BOOL "" FORCE)
tribits_project()

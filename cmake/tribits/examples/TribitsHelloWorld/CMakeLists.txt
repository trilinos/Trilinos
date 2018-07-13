# Example of a simple project that uses TriBITS
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11 FATAL_ERROR)
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/ProjectName.cmake)
PROJECT(${PROJECT_NAME} NONE)
SET(${PROJECT_NAME}_TRIBITS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../..
  CACHE PATH "TriBITS base directory (default assumes in TriBITS source tree).")
INCLUDE("${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake")
# Only one package in this simple project so just enable it :-)
SET(${PROJECT_NAME}_ENABLE_HelloWorld ON CACHE BOOL "" FORCE)
TRIBITS_PROJECT()

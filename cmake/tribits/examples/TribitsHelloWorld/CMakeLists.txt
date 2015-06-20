# Example of a simple project that uses TriBITS

# To be safe, define your minimum CMake version.  This may be newer than the
# min required by TriBITS.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11 FATAL_ERROR)

# Make CMake set WIN32 with CYGWIN for older CMake versions.  CMake requires
# this to be in the top-level CMakeLists.txt file and not an include file :-(
SET(CMAKE_LEGACY_CYGWIN_WIN32 1 CACHE BOOL "" FORCE)

# Get PROJECT_NAME (must be in file for other parts of system)
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/ProjectName.cmake)

# CMake requires that you declare the CMake project in the top-level file and
# not in an include file :-(
PROJECT(${PROJECT_NAME} NONE)

SET(${PROJECT_NAME}_TRIBITS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../..
  CACHE PATH "TriBITS base directory (default assumes in TriBITS source tree).")

# Include the TriBITS system
INCLUDE("${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake")

# Disable MPI and Fortran
SET(TPL_ENABLE_MPI OFF CACHE BOOL "" FORCE)
# Turn off Fortran support by default
SET(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT OFF)
# Turn of export makefiles and the like
SET(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES_DEFAULT OFF)
# Only one package in this simple project so just enable it :-)
SET(${PROJECT_NAME}_ENABLE_HelloWorld ON CACHE BOOL "" FORCE)

# Do all of the processing for this Tribits project
TRIBITS_PROJECT()

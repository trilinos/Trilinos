#################################################################################
#                                                                               #
# SEACAS -- Build as package in Trilinos or as standalone project using Tribits # 
#                                                                               #
#################################################################################

macro(SEACAS_PACKAGE target)
  if(${CMAKE_PROJECT_NAME} STREQUAL "SEACAS")
    TRIBITS_PACKAGE(${target})
  else()
    INCLUDE(TribitsSubPackageMacros)
    TRIBITS_SUBPACKAGE(${target})
  endif()
endmacro()

macro(SEACAS_POSTPROCESS)
  if(${CMAKE_PROJECT_NAME} STREQUAL "SEACAS")
    TRIBITS_PACKAGE_POSTPROCESS()
  else()
    TRIBITS_SUBPACKAGE_POSTPROCESS()
  endif()
endmacro()

function(install_executable target)
  INSTALL(TARGETS ${target} COMPONENT ${PACKAGE_NAME} RUNTIME DESTINATION ${${PROJECT_NAME}_INSTALL_RUNTIME_DIR})
endfunction()

if(NOT CMAKE_PROJECT_NAME OR "${CMAKE_PROJECT_NAME}" STREQUAL "SEACAS")
  #
  # This is being built as a standalone project
  #
  MESSAGE("----Building SEACAS standalone")
  cmake_minimum_required(VERSION 3.0)
  # Make CMake set WIN32 with CYGWIN for older CMake versions.  CMake requires
  # this to be in the top-level CMakeLists.txt file and not an include file :-(
  SET(CMAKE_LEGACY_CYGWIN_WIN32 1 CACHE BOOL "" FORCE)

  # A) Define your project name and set up major project options
  # Get PROJECT_NAME (must be in a file for other parts of system to read)
  INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/ProjectName.cmake")

  # CMake requires that you declare the CMake project in the top-level file and
  # not in an include file :-(
  PROJECT(${PROJECT_NAME} NONE)

  SET(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT TRUE)

  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cmake/tribits/TriBITS.cmake)
  SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

  SET(SUBDIR "src/")

  # Used for dependencies SEACASExodus vs Exodus
  SET(SEACAS_PREFIX "")  
  SET(IOSS_PREFIX "IOSS")
else()
  #
  # This is being built as a package in Trilinos
  #
  MESSAGE("----Building SEACAS as package in ${CMAKE_PROJECT_NAME}")

  TRIBITS_PACKAGE_DECL( SEACAS ENABLE_SHADOWING_WARNINGS DISABLE_STRONG_WARNINGS )

  TRIBITS_ADD_DEBUG_OPTION()
  TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

  TRIBITS_PROCESS_SUBPACKAGES()

  TRIBITS_PACKAGE_DEF()
  TRIBITS_PACKAGE_POSTPROCESS()

  TRIBITS_EXCLUDE_FILES(doc)

  SET(SUBDIR "")

  # Used for dependencies SEACASExodus vs Exodus
  SET(SEACAS_PREFIX "SEACAS")
  SET(IOSS_PREFIX "SEACAS")
endif()

INCLUDE(TribitsGlobalMacros)
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

if(${CMAKE_PROJECT_NAME} STREQUAL "SEACAS")
  # Do all of the processing for this Tribits project
  TRIBITS_PROJECT()
  # If don't clear PACKAGE_NAME, it ends up being set to last entry in PackagesList.cmake
  SET(PACKAGE_NAME "")  
endif()

# For decomp, blot, fastq, and exodus.py
ADD_SUBDIRECTORY(${SUBDIR}scripts)



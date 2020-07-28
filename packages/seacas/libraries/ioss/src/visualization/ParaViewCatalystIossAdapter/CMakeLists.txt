# Copyright(C) 1999-2020 National Technology & Engineering Solutions
# of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
# NTESS, the U.S. Government retains certain rights in this software.
#
# See packages/seacas/LICENSE for details

cmake_minimum_required (VERSION 2.8.1)
project(ParaViewCatalystIossAdapter)

set(PVCSA_VERSION CACHE STRING "Version string to append to catalyst adapter library.")

if(PVCSA_VERSION)
 SET(PVCSA_LIBRARY_NAME ParaViewCatalystIossAdapter_${PVCSA_VERSION})
else()
 SET(PVCSA_LIBRARY_NAME ParaViewCatalystIossAdapter)
endif()

find_package(ParaView REQUIRED COMPONENTS vtkPVCatalyst vtkPVPythonCatalyst)
include(${PARAVIEW_USE_FILE})

set(PARAVIEW_PYTHON_LIBRARY_RPATH_DIRECTORY CACHE
    FILEPATH "Path to the directory containing the Python library used by ParaView.")

if(PARAVIEW_PYTHON_LIBRARY_RPATH_DIRECTORY)
  cmake_policy(SET CMP0015 OLD)
  link_directories(${PARAVIEW_PYTHON_LIBRARY_RPATH_DIRECTORY})
endif()

SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" ${PARAVIEW_PYTHON_LIBRARY_RPATH_DIRECTORY})
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

include_directories("include")

add_subdirectory(parser)
add_library(${PVCSA_LIBRARY_NAME} SHARED ParaViewCatalystIossAdapter.cxx
                                         vtkExodusIIMultiBlockDataSet.cxx)
target_link_libraries(${PVCSA_LIBRARY_NAME} vtkPVCatalyst)
target_link_libraries(${PVCSA_LIBRARY_NAME} vtkPVPythonCatalyst)
target_link_libraries(${PVCSA_LIBRARY_NAME} catalystparser)

INSTALL(TARGETS ${PVCSA_LIBRARY_NAME} LIBRARY DESTINATION lib)
INSTALL(DIRECTORY python DESTINATION .)

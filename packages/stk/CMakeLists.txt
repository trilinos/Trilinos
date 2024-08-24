# Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
# Solutions of Sandia, LLC (NTESS). Under the terms of Contract
# DE-NA0003525 with NTESS, the U.S. Government retains certain rights
# in this software.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
# 
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
# 
#     * Neither the name of NTESS nor the names of its contributors
#       may be used to endorse or promote products derived from this
#       software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 

cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

message("starting STK cmake configuration, CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}")

option(STK_BUILT_FOR_SIERRA "Enable SIERRA capability" OFF)
set(SIERRA_MIGRATION ${STK_BUILT_FOR_SIERRA} CACHE BOOL "Enable SIERRA capability")

if(STK_BUILT_FOR_SIERRA)
set(STK_ENABLE_16BIT_UPWARDCONN_INDEX_TYPE ON)
endif()

IF(COMMAND TRIBITS_PACKAGE_DECL)
  SET(HAVE_STK_Trilinos ON)
  TRIBITS_PACKAGE_DECL(STK)
  MESSAGE("*** Building STK as a Trilinos package. ***")
ELSE()
  SET(HAVE_STK_Trilinos OFF)
  project(STK CXX Fortran)

  SET(PACKAGE_NAME "STK")

  MESSAGE("*** Building STK as a stand-alone cmake package. ***")
ENDIF()
  
SET(STK_TOPLEVEL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(STK_TOPLEVEL_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

INCLUDE(${STK_TOPLEVEL_SOURCE_DIR}/cmake/stk_wrappers.cmake)

if(CMAKE_INSTALL_INCLUDE_DIR)
  SET(STK_INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDE_DIR})
elseif(INSTALL_INCLUDE_DIR)
  SET(STK_INSTALL_INCLUDE_DIR ${INSTALL_INCLUDE_DIR})
else()
  SET(STK_INSTALL_INCLUDE_DIR "include")
endif()

if(CMAKE_INSTALL_LIBDIR)
  SET(STK_INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR})
elseif(INSTALL_LIBDIR)
  SET(STK_INSTALL_LIBDIR ${INSTALL_LIBDIR})
else()
  SET(STK_INSTALL_LIBDIR "lib")
endif()

if(CMAKE_INSTALL_BINDIR)
  SET(STK_INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR})
elseif(INSTALL_BINDIR)
  SET(STK_INSTALL_BINDIR ${INSTALL_BINDIR})
else()
  SET(STK_INSTALL_BINDIR "bin")
endif()

STK_ADD_DEBUG_AND_DEPRECATED_OPTIONS()

MESSAGE("\nPROJECT_NAME: '${PROJECT_NAME}'")
MESSAGE("PACKAGE_NAME: '${PACKAGE_NAME}'")
MESSAGE("${PACKAGE_NAME}_SOURCE_DIR: '${${PACKAGE_NAME}_SOURCE_DIR}'\n")

IF (HAVE_STK_Trilinos)
  SET(STK_HAVE_KOKKOS ON)
  if (TPL_ENABLE_MPI OR HAVE_MPI)
    message("STK_HAS_MPI is true")
    SET(STK_HAS_MPI ON)
  else()
    message("STK_HAS_MPI is false")
    SET(STK_HAS_MPI OFF)
  endif()
ELSE()
  if(STK_ENABLE_TESTS)
    enable_testing()
  endif()

  find_package(Kokkos REQUIRED)
  SET(STK_HAVE_KOKKOS ON)

  IF(DEFINED STK_ENABLE_MPI)
    IF(STK_ENABLE_MPI)
      MESSAGE("MPI requested via STK_ENABLE_MPI=ON")
    ELSE()
      MESSAGE("MPI disabled via STK_ENABLE_MPI=OFF")
    ENDIF()
  ELSE()
    MESSAGE("MPI defaulting to off. (STK_ENABLE_MPI not defined)")
  ENDIF()

  IF(STK_ENABLE_MPI)
    find_package(MPI)
    IF(MPI_FOUND)
      SET(STK_HAS_MPI ON)
      include_directories(SYSTEM ${MPI_INCLUDE_PATH})
      message("MPI_INCLUDE_PATH: ${MPI_INCLUDE_PATH}")
    ELSE()
      MESSAGE(FATAL_ERROR "MPI enabled by '-DSTK_ENABLE_MPI' but not found.")
    ENDIF()
  ELSE()
    MESSAGE("Building serial without MPI. (To enable MPI, use '-DSTK_ENABLE_MPI:BOOL=ON')")
  ENDIF()
ENDIF()

find_package(ArborX QUIET)
if(TARGET ArborX::ArborX)
  MESSAGE("Found ArborX, making it available within stk")
  SET(STK_HAS_ARBORX ON)
else()
  MESSAGE("Optional search library ArborX is not enabled.")
endif()

if(NOT HAVE_STK_Trilinos)
  find_package(SEACAS)
endif()

if(SEACAS_ENABLE_SEACASIoss)
  message("Enabling stk usage of SEACASIoss")
  SET(STK_HAS_SEACAS_IOSS ON)
endif()
if(SEACAS_ENABLE_SEACASExodus)
  message("Enabling stk usage of SEACASExodus")
  SET(STK_HAS_SEACAS_EXODUS ON)
endif()
if(SEACAS_ENABLE_SEACASNemesis)
  message("Enabling stk usage of SEACASNemesis")
  SET(STK_HAS_SEACAS_NEMESIS ON)
endif()

if (Trilinos_ENABLE_Intrepid2)
  message("Intrepid2 is enabled, setting STK_HAVE_INTREPID2")
  SET(STK_HAVE_INTREPID2 ON)
endif()

if(NOT HAVE_STK_Trilinos)
  stk_process_enables()
endif()

if(STK_ENABLE_16BIT_UPWARDCONN_INDEX_TYPE)
  set(STK_16BIT_UPWARDCONN_INDEX_TYPE ON)
  message("STK_ENABLE_16BIT_UPWARDCONN_INDEX_TYPE is ON, upward-connectivity tables limited to ~65K entries per bucket.")
else()
  set(STK_16BIT_UPWARDCONN_INDEX_TYPE OFF)
  message("STK_ENABLE_16BIT_UPWARDCONN_INDEX_TYPE is OFF, upward-connectivity tables using 32bit indexing, (allows ~4B entries per bucket).")
endif()

IF (TPL_ENABLE_Boost)
  SET(STK_HAVE_BOOST ON)
ENDIF()
IF (${${PROJECT_NAME}_ENABLE_STKNGP_TEST})
  SET(STK_HAVE_STKNGP_TEST ON)
ENDIF()
IF (${${PROJECT_NAME}_ENABLE_STKMesh})
  SET(STK_HAVE_STKMESH ON)
ENDIF()
IF (${${PROJECT_NAME}_ENABLE_STKIO})
  SET(STK_HAVE_STKIO ON)
ENDIF()

SET(STK_VOLATILE_SIMD ON)

IF ( ${PACKAGE_NAME}_DISABLE_MPI_NEIGHBOR_COMM )
  SET(STK_DISABLE_MPI_NEIGHBOR_COMM TRUE)
ENDIF()

IF (NOT ${PROJECT_NAME}_ENABLE_SEACASAprepro_lib)
  SET(NOT_HAVE_STK_SEACASAPREPRO_LIB TRUE)
ENDIF()

# set fortran mangling options

IF("${FC_FN_UNDERSCORE}" STREQUAL "NO_UNDER")
  SET(FORTRAN_NO_UNDERSCORE ON)
ELSEIF("${FC_FN_UNDERSCORE}" STREQUAL "UNDER")
  SET(FORTRAN_ONE_UNDERSCORE ON)
ELSEIF("${FC_FN_UNDERSCORE}" STREQUAL "SECOND_UNDER")
  SET(FORTRAN_TWO_UNDERSCORES ON)
ELSE()
  MESSAGE("Could not determine the Fortran mangling; defaulting to one underscore.")
  SET(FORTRAN_ONE_UNDERSCORE ON)
ENDIF()

MESSAGE("\nCMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

STK_CONFIGURE_FILE(STK_Trilinos_config.h)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/STK_Trilinos_config.h DESTINATION
        ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_util)

STK_SUBPACKAGES()

STK_PACKAGE_POSTPROCESS()

IF (NOT HAVE_STK_Trilinos)
  include(CMakePackageConfigHelpers)
  
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/stkConfigVersion.cmake"
    VERSION VOTD
    COMPATIBILITY ExactVersion
  )

  install(EXPORT stkTargets
    FILE stkTargets.cmake
    NAMESPACE stk::
    DESTINATION share/cmake/stk
  )

  install(
    FILES
      cmake/stkConfig.cmake
      "${CMAKE_CURRENT_BINARY_DIR}/stkConfigVersion.cmake"
    DESTINATION share/cmake/stk
    COMPONENT Devel
  )
  if(STK_BUILT_FOR_SIERRA)
    install(
      FILES cmake/stkLapackSierra.cmake
      DESTINATION share/cmake/stk
      RENAME stkLapack.cmake
      COMPONENT Devel
    )
  else()
    install(
      FILES cmake/stkLapackGeneric.cmake
      DESTINATION share/cmake/stk
      RENAME stkLapack.cmake
      COMPONENT Devel
    )
  endif()
ENDIF()

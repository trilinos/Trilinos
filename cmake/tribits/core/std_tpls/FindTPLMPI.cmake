# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

# Either the MPI compiler wrappers take care of these or the user has to set
# the explicitly using basic compile flags and ${PROJECT_NAME}_EXTRA_LINK_FLAGS.
global_set(TPL_MPI_INCLUDE_DIRS "")
global_set(TPL_MPI_LIBRARIES "")
global_set(TPL_MPI_LIBRARY_DIRS "")

if(WIN32 AND TPL_ENABLE_MPI)
  find_package(MPI)
  include_directories(${MPI_INCLUDE_PATH})
  global_set(TPL_MPI_INCLUDE_DIRS ${MPI_INCLUDE_PATH})
  # FindMPI often leaves MPI_LIBRARIES empty on Windows while MPI_C/CXX_LIBRARIES
  # hold msmpi.lib.  Also, global_set() does not use FORCE, so avoid leaving a
  # stale empty TPL_MPI_LIBRARIES in the cache from an earlier configure.
  set(_tpl_mpi_link_libs "")
  if(MPI_LIBRARIES)
    set(_tpl_mpi_link_libs "${MPI_LIBRARIES}")
  elseif(MPI_CXX_LIBRARIES)
    set(_tpl_mpi_link_libs "${MPI_CXX_LIBRARIES}")
  elseif(MPI_C_LIBRARIES)
    set(_tpl_mpi_link_libs "${MPI_C_LIBRARIES}")
  endif()
  if(_tpl_mpi_link_libs)
    set(TPL_MPI_LIBRARIES "${_tpl_mpi_link_libs}" CACHE INTERNAL "" FORCE)
  endif()
  # When MPI_LIBRARY_NAMES is empty, tribits_tpl_find_include_dirs_and_libraries(MPI)
  # calls global_null_set(TPL_MPI_LIBRARIES) and clears the cache.  If FindMPI
  # already produced MPI::MPI_C / MPI::MPI_CXX, create MPI::all_libs here so
  # that generic finder is skipped.
  if(MPI_C_FOUND AND MPI_CXX_FOUND AND NOT TARGET MPI::all_libs)
    tribits_extpkg_create_imported_all_libs_target_and_config_file(
      MPI
      INNER_FIND_PACKAGE_NAME MPI
      IMPORTED_TARGETS_FOR_ALL_LIBS MPI::MPI_C MPI::MPI_CXX)
  endif()
endif()

# Don't allow calling find_package(MPI) by default.  Force the user to set
# MPI_ALLOW_PACKAGE_PREFIND=TRUE or even MPI_FORCE_PRE_FIND_PACKAGE if then
# want to force the finding of MPI using find_package(MPI).
set(MPI_ALLOW_PACKAGE_PREFIND  FALSE  CACHE  BOOL
  "Allow calling find_package(MPI) by default (default is FALSE)")

tribits_tpl_allow_pre_find_package(MPI MPI_ALLOW_PREFIND)
if(MPI_ALLOW_PREFIND)
  find_package(MPI)
  if(MPI_C_FOUND AND MPI_CXX_FOUND)
    tribits_extpkg_create_imported_all_libs_target_and_config_file(
      MPI
      INNER_FIND_PACKAGE_NAME MPI
      IMPORTED_TARGETS_FOR_ALL_LIBS MPI::MPI_C MPI::MPI_CXX)
  endif()
endif()

if(NOT TARGET MPI::all_libs)
  tribits_tpl_find_include_dirs_and_libraries(MPI)
  # NOTE: Above, we need to generate the MPI::all_libs target and the
  # MPIConfig.cmake file that will also provide the MPI::all_libs target.
endif()

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
  global_set(TPL_MPI_LIBRARIES ${MPI_LIBRARIES})
endif()

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
endif()

# NOTE: Above, we need to generate the MPI::all_libs target and the
# MPIConfig.cmake file that will also provide the MPI::all_libs target.

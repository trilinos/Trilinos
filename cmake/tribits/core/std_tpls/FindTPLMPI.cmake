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

tribits_tpl_find_include_dirs_and_libraries(MPI)

# NOTE: Above, we need to generate the MPI::all_libs target and the
# MPIConfig.cmake file that will also provide the MPI::all_libs target.

# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

get_filename_component(_Trilinos_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(_Trilinos_PREFIX "${_Trilinos_DIR}" PATH)
set(Trilinos_DIR "${_Trilinos_PREFIX}/lib/cmake/Trilinos")
message(WARNING "TrilinosConfig.cmake has moved.  "
  "It now exists at a location under the installation prefix where the "
  "find_package command looks by default (<prefix>/lib/cmake/Trilinos).  "
  "This compatibility file exists at the old location (<prefix>/include) "
  "to present this message and load the file from its new location."
  "\n"
  "The find_package() call that loaded this file did so because its "
  "cached result variable, Trilinos_DIR, is set to\n"
  "  ${_Trilinos_DIR}\n"
  "I'm locally setting Trilinos_DIR to\n"
  "  ${Trilinos_DIR}\n"
  "and loading TrilinosConfig.cmake from its new location.  "
  "One may suppress this warning by setting the above value in the cache.  "
  "However, the application needs modification permanently fix the issue.  "
  "The find_package() call that loaded this file may have the form\n"
  "  find_package(Trilinos REQUIRED PATHS \${Trilinos_PATH}/include)\n"
  "Change it to the form\n"
  "  set(CMAKE_PREFIX_PATH \${Trilinos_PATH} \${CMAKE_PREFIX_PATH})\n"
  "  find_package(Trilinos REQUIRED)\n"
  "to find TrilinosConfig.cmake in its new location in future builds "
  "while still honoring the Trilinos_PATH option for this application."
  )
include(${Trilinos_DIR}/TrilinosConfig.cmake)

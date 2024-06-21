# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

#
# First, set up the variables for the (backward-compatible) TriBITS way of
# finding Netcdf.  These are used in case find_package(NetCDF ...) is not
# called or does not find NetCDF.  Also, these variables need to be non-null
# in order to trigger the right behavior in the function
# tribits_tpl_find_include_dirs_and_libraries().
#
set(REQUIRED_HEADERS cgnslib.h)
set(REQUIRED_LIBS_NAMES cgns)

#
# Second, search for Netcdf components (if allowed) using the standard
# find_package(CGNS ...).
#
tribits_tpl_allow_pre_find_package(CGNS  CGNS_ALLOW_PREFIND)
if (CGNS_ALLOW_PREFIND)

  message("-- Using find_package(CGNS ...) ...")

  set(CMAKE_MODULE_PATH
    "${CMAKE_MODULE_PATH}"
    "${CMAKE_CURRENT_LIST_DIR}/find_modules"
    "${CMAKE_CURRENT_LIST_DIR}/utils"
     )
  
  find_package(CGNS)

  if (CGNS_FOUND)
    set(TPL_CGNS_LIBRARIES ${CGNS_LIBRARIES} CACHE PATH
      "List of semi-colon separated (full) paths to the CGNS libraries")
    set(TPL_CGNS_INCLUDE_DIRS ${CGNS_INCLUDE_DIRS} CACHE PATH
      "List of semi-colon separated list of directories containing CGNS header files")
  endif()

endif()

#
# Third, call tribits_tpl_find_include_dirs_and_libraries()
#
tribits_tpl_find_include_dirs_and_libraries( CGNS
  REQUIRED_HEADERS ${REQUIRED_HEADERS}
  REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES})

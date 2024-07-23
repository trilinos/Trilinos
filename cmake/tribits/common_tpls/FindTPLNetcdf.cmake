# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

cmake_policy(SET CMP0074 NEW)

set(Netcdf_ALLOW_MODERN FALSE CACHE BOOL "Allow finding Netcdf as a modern CMake config file with exported targets (and only this way)")

if (Netcdf_ALLOW_MODERN AND HDF5_FOUND_MODERN_CONFIG_FILE)

  set(minimum_modern_netCDF_version 4.7.4)
  print_var(Netcdf_ALLOW_MODERN)
  message("-- Using find_package(netCDF ${minimum_modern_netCDF_version} CONFIG) ...")
  find_package(netCDF ${minimum_modern_netCDF_version} CONFIG)
  if (netCDF_FOUND)
    message("-- Found netCDF_CONFIG=${netCDF_CONFIG}")
    message("-- Generating Netcdf::all_libs and NetcdfConfig.cmake")
    tribits_extpkg_create_imported_all_libs_target_and_config_file(
      Netcdf
      INNER_FIND_PACKAGE_NAME  netCDF
      IMPORTED_TARGETS_FOR_ALL_LIBS   netCDF::netcdf)
    set(TPL_Netcdf_NOT_FOUND FALSE)
  endif()

endif()

if (NOT TARGET  Netcdf::all_libs)

  # First, set up the variables for the (backward-compatible) TriBITS way of
  # finding Netcdf.  These are used in case find_package(Netcdf ...) is not
  # called or does not find HDF5.  Also, these variables need to be non-null
  # in order to trigger the right behavior in the function
  # tribits_tpl_find_include_dirs_and_libraries().

  set(REQUIRED_HEADERS netcdf.h)
  set(REQUIRED_LIBS_NAMES netcdf)

  if (TPL_ENABLE_MPI)
    set(REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES} pnetcdf)
  endif()

  #
  # Second, search for Netcdf components (if allowed) using the standard
  # find_package(NetCDF ...).
  #
  tribits_tpl_allow_pre_find_package(Netcdf  Netcdf_ALLOW_PREFIND)
  if (Netcdf_ALLOW_PREFIND)

    message("-- Using find_package(Netcdf ...) ...")

    set(CMAKE_MODULE_PATH
      "${CMAKE_MODULE_PATH}"
      "${CMAKE_CURRENT_LIST_DIR}/find_modules"
      "${CMAKE_CURRENT_LIST_DIR}/utils"
       )

    find_package(NetCDF)

    if (NetCDF_FOUND)
      set(DOCSTR "List of semi-colon separated paths to look for the TPL Netcdf")
      set(TPL_Netcdf_Enables_Netcdf4 ${NetCDF_NEEDS_HDF5} CACHE BOOL
        "True if netcdf enables netcdf-4")
      set(TPL_Netcdf_Enables_PNetcdf ${NetCDF_NEEDS_PNetCDF} CACHE BOOL
        "True if netcdf enables pnetcdf")
      set(TPL_Netcdf_PARALLEL ${NetCDF_PARALLEL} CACHE INTERNAL
        "True if netcdf compiled with parallel enabled")
      set(TPL_Netcdf_LIBRARY_DIRS ${_hdf5_LIBRARY_SEARCH_DIRS} CACHE PATH
        "${DOCSTR} library files")
      set(TPL_Netcdf_LIBRARIES ${NetCDF_LIBRARIES} CACHE PATH
        "List of semi-colon separated library names (not 'lib' or extension).")
      set(TPL_Netcdf_INCLUDE_DIRS ${NetCDF_INCLUDE_DIRS} CACHE PATH
        "${DOCSTR} header files.")
    endif()
  else()
    # Curl library is only required if DAP is enabled; should detect inside
    # FindNetCDF.cmake, but that is not being called... SEMS has DAP enabled;
    # many HPC systems don't, but they override these settings...
    find_program(NC_CONFIG "nc-config")
    if (NC_CONFIG)
      execute_process(COMMAND "nc-config --has-dap2"
                      OUTPUT_VARIABLE NETCDF_HAS_DAP2)
      execute_process(COMMAND "nc-config --has-dap4"
                      OUTPUT_VARIABLE NETCDF_HAS_DAP4)
    endif()
    if ((NOT NC_CONFIG) OR NETCDF_HAS_DAP2 OR NETCDF_HAS_DAP4)
      set(REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES} curl)
    endif()
  endif()

  #
  # Third, call tribits_tpl_find_include_dirs_and_libraries()
  #
  tribits_tpl_find_include_dirs_and_libraries( Netcdf
    REQUIRED_HEADERS ${REQUIRED_HEADERS}
    REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES}
    )
  # NOTE: If find_package(Netcdf ...) was called and successfully found Netcdf,
  # then tribits_tpl_find_include_dirs_and_libraries() will use the already-set
  # variables TPL_Netcdf_INCLUDE_DIRS and TPL_Netcdf_LIBRARIES and then print
  # them out (and set some other standard variables as well).  This is the final
  # "hook" into the TriBITS TPL system.

  # If the `find_package(NetCDF)` is not run, then this may not be set
  # Need to determine how this is set in the library that is being used...

  if ("${TPL_Netcdf_PARALLEL}" STREQUAL "")
     assert_defined(TPL_Netcdf_INCLUDE_DIRS)
     find_path(meta_path
        NAMES "netcdf_meta.h"
        HINTS ${TPL_Netcdf_INCLUDE_DIRS}
        NO_DEFAULT_PATH)

     if (meta_path)
        # Search meta for NC_HAS_PARALLEL setting...
        # Note that there is both NC_HAS_PARALLEL and NC_HAS_PARALLEL4, only want first...
        file(STRINGS "${meta_path}/netcdf_meta.h" netcdf_par_string REGEX "NC_HAS_PARALLEL ")
        string(REGEX MATCH "[01]" netcdf_par_val "${netcdf_par_string}")
        if (netcdf_par_val EQUAL 1)
           set(TPL_Netcdf_PARALLEL True CACHE INTERNAL
               "True if netcdf compiled with parallel enabled")
        endif()
     endif()
     if ("${TPL_Netcdf_PARALLEL}" STREQUAL "")
        set(TPL_Netcdf_PARALLEL False CACHE INTERNAL
            "True if netcdf compiled with parallel enabled")
     endif()
  endif()
  message(STATUS "TPL_Netcdf_PARALLEL is ${TPL_Netcdf_PARALLEL}")
endif()

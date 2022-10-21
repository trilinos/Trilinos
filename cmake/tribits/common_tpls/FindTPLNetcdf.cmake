# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2016 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER

#
# First, set up the variables for the (backward-compatible) TriBITS way of
# finding Netcdf.  These are used in case find_package(NetCDF ...) is not
# called or does not find NetCDF.  Also, these variables need to be non-null
# in order to trigger the right behavior in the function
# tribits_tpl_find_include_dirs_and_libraries().
#
if (${CMAKE_VERSION} GREATER "3.13")
     cmake_policy(SET CMP0074 NEW)
endif()

set(Netcdf_ALLOW_MODERN FALSE CACHE BOOL "Allow finding Netcdf as a modern CMake config file with exported targets (and only this way)")

if (Netcdf_ALLOW_MODERN)

  set(minimum_modern_netCDF_version 4.7.4)
  message("-- Netcdf_ALLOW_MODERN=${Netcdf_ALLOW_MODERN}")
  message("-- Using find_package(netCDF ${minimum_modern_netCDF_version} CONFIG) ...")
  find_package(netCDF ${minimum_modern_netCDF_version} CONFIG)
  if (netCDF_FOUND)
    message("-- Found netCDF_CONFIG=${netCDF_CONFIG}")
    message("-- Generating Netcdf::all_libs and NetcdfConfig.cmake")
    # instead of using tribits_extpkg_create_imported_all_libs_target_and_config_file,
    # we will do what it does ourselves because we need to bring in two
    # inner find packages (netCDF and hdf5) because netCDF does not properly
    # do find_dependency(hdf5)
    add_library(Netcdf::all_libs  INTERFACE  IMPORTED  GLOBAL)
    target_link_libraries(Netcdf::all_libs  INTERFACE  netCDF::netcdf)
    set(configFileStr "")
    string(APPEND configFileStr
      "include(CMakeFindDependencyMacro)\n" )
    string(APPEND configFileStr
      "set(netCDF_DIR \"${netCDF_DIR}\")\n" )
    string(APPEND configFileStr
      "find_dependency(netCDF ${minimum_modern_netCDF_version} CONFIG)\n"
      )
    if (netCDF_HAS_HDF5)
      find_package(hdf5 CONFIG REQUIRED)
      string(APPEND configFileStr
        "set(hdf5_DIR \"${hdf5_DIR}\")\n" )
      string(APPEND configFileStr
        "find_dependency(hdf5 CONFIG)\n"
        )
    endif()
    string(APPEND configFileStr
      "add_library(Netcdf::all_libs  INTERFACE  IMPORTED  GLOBAL)\n"
      )
    string(APPEND configFileStr
      "target_link_libraries(Netcdf::all_libs  INTERFACE  netCDF::netcdf)\n")
    set(buildDirExternalPkgsDir
      "${${PROJECT_NAME}_BINARY_DIR}/${${PROJECT_NAME}_BUILD_DIR_EXTERNAL_PKGS_DIR}")
    set(tplConfigFile
      "${buildDirExternalPkgsDir}/Netcdf/NetcdfConfig.cmake")
    file(WRITE "${tplConfigFile}" "${configFileStr}")
    # the following gets set by tribits_tpl_find_include_dirs_and_libraries
    # and it is not well-documented that you need to set it yourself if not using
    # tribits_tpl_find_include_dirs_and_libraries
    set(TPL_Netcdf_NOT_FOUND FALSE)
    message("-- TPL_Netcdf_PARALLEL = netCDF_HAS_PARALLEL = ${netCDF_HAS_PARALLEL}")
    set(TPL_Netcdf_PARALLEL ${netCDF_HAS_PARALLEL})
    message("-- TPL_Netcdf_Enables_Netcdf4 = netCDF_HAS_NC4 = ${netCDF_HAS_NC4}")
    set(TPL_Netcdf_Enables_Netcdf4 ${netCDF_HAS_NC4} CACHE BOOL
      "True if netcdf enables netcdf-4")
  endif()

else()

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

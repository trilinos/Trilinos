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
# finding Netcdf.  These are used in case FIND_PACKAGE(NetCDF ...) is not
# called or does not find NetCDF.  Also, these variables need to be non-null
# in order to trigger the right behavior in the function
# TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES().
#

SET(REQUIRED_HEADERS netcdf.h)
SET(REQUIRED_LIBS_NAMES netcdf)

#
# Second, search for Netcdf components (if allowed) using the standard
# FIND_PACKAGE(NetCDF ...).
#
TRIBITS_TPL_ALLOW_PRE_FIND_PACKAGE(Netcdf  Netcdf_ALLOW_PREFIND)
IF (Netcdf_ALLOW_PREFIND)

  MESSAGE("-- Using FIND_PACKAGE(Netcdf ...) ...")

  SET(CMAKE_MODULE_PATH
    "${CMAKE_MODULE_PATH}"
    "${CMAKE_CURRENT_LIST_DIR}/find_modules"
    "${CMAKE_CURRENT_LIST_DIR}/utils"
     )
  
  find_package(NetCDF)

  IF (NetCDF_FOUND)
    SET(DOCSTR "List of semi-colon separated paths to look for the TPL Netcdf")
    set(TPL_Netcdf_Enables_Netcdf4 ${NetCDF_NEEDS_HDF5} CACHE BOOL
      "True if netcdf enables netcdf-4")
    set(TPL_Netcdf_Enables_PNetcdf ${NetCDF_NEEDS_PNetCDF} CACHE BOOL
      "True if netcdf enables pnetcdf")
    set(TPL_Netcdf_PARALLEL ${NetCDF_PARALLEL} CACHE BOOL
      "True if netcdf compiled with parallel enabled")
    set(TPL_Netcdf_LIBRARY_DIRS ${_hdf5_LIBRARY_SEARCH_DIRS} CACHE PATH
      "${DOCSTR} library files")
    set(TPL_Netcdf_LIBRARIES ${NetCDF_LIBRARIES} CACHE PATH
      "List of semi-colon seprated library names (not 'lib' or extension).")
    set(TPL_Netcdf_INCLUDE_DIRS ${NetCDF_INCLUDE_DIRS} CACHE PATH
      "${DOCSTR} header files.")
  ENDIF()

ENDIF()

#
# Third, call TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES()
#
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( Netcdf
  REQUIRED_HEADERS ${REQUIRED_HEADERS}
  REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES}
  )
# NOTE: If FIND_PACKAGE(Netcdf ...) was called and successfully found Netcdf,
# then TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES() will use the already-set
# variables TPL_Netcdf_INCLUDE_DIRS and TPL_Netcdf_LIBRARIES and then print
# them out (and set some other standard variables as well).  This is the final
# "hook" into the TriBITS TPL system.

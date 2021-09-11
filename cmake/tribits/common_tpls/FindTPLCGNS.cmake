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

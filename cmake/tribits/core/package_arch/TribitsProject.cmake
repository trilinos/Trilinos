# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
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

# Projects that change the location of the source need to consider this in
# their top-level CMakeLists.txt file
set(${PROJECT_NAME}_TRIBITS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/tribits"
  CACHE PATH
  "The base directory pointing to the TriBITS system.  If provided as a relative path (formatted as STRING) then will be set to '${CMAKE_CURRENT_SOURCE_DIR}/'.  NOTE: If you leave off the STRING datatype, and it is a relative path, then it will be interpreted as relative to the build directory!"
  )
mark_as_advanced(${PROJECT_NAME}_TRIBITS_DIR)

if (NOT IS_ABSOLUTE "${${PROJECT_NAME}_TRIBITS_DIR}")
  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("NOTE: ${PROJECT_NAME}_TRIBITS_DIR = '${${PROJECT_NAME}_TRIBITS_DIR}' provided as a relative directory so making this relative to '${CMAKE_CURRENT_SOURCE_DIR}'!")
  endif()
  set(${PROJECT_NAME}_TRIBITS_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/${${PROJECT_NAME}_TRIBITS_DIR}")
endif()

if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
  message("${PROJECT_NAME}_TRIBITS_DIR='${${PROJECT_NAME}_TRIBITS_DIR}'")
endif()

set(CMAKE_MODULE_PATH
   ${${PROJECT_NAME}_TRIBITS_DIR}/core/package_arch
   )

if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
  message("CMAKE_MODULE_PATH='${CMAKE_MODULE_PATH}'")
endif()

# Overrides that we have for CMake functions
include(TribitsCMakePolicies  NO_POLICY_SCOPE)
include(TribitsProjectImpl)


# @MACRO: tribits_project()
#
# Processes a `TriBITS Project`_'s files and configures its software which is
# called from the project's top-level `<projectDir>/CMakeLists.txt`_ file.
#
# Usage::
#
#   tribits_project()
#
# This macro requires that the variable `PROJECT_NAME`_ be defined before
# calling this macro.  All default values for project settings should be set
# before calling this macro (see `TriBITS Global Project Settings`_).  Also,
# the variable `${PROJECT_NAME}_TRIBITS_DIR`_ must be set as well.
#
# This macro then adds all of the necessary paths to ``CMAKE_MODULE_PATH`` and
# then performs all processing of the TriBITS project files (see `Full TriBITS
# Project Configuration`_).
#
macro(tribits_project)
  tribits_project_impl(${ARGN})
endmacro()

# Note, this is just a shell of a macro that calls the real implementation
# tribits_project_impl().  This allows someone to set
# ${PROJECT_NAME}_TRIBITS_DIR in the env and point to a different Tribits
# implementation to test before snapshoting.

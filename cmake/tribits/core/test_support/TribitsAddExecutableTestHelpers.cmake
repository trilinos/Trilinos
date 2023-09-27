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

include_guard()

include("${CMAKE_CURRENT_LIST_DIR}/../common/TribitsCMakePolicies.cmake"  NO_POLICY_SCOPE)

include("${CMAKE_CURRENT_LIST_DIR}/../utils/AdvancedSet.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/../utils/MessageWrapper.cmake")

advanced_set( ${PROJECT_NAME}_CMAKE_EXECUTABLE_SUFFIX ".exe"
  CACHE STRING
  "Default exec suffix on all platforms (can be overridden by each executable added)." )


# Process the COMM arguments
#
# NOTE: The COMM array arguments is passed as ${ARGN}
#
function(tribits_process_comm_args  ADD_SERIAL_FEATURE_OUT  ADD_MPI_FEATURE_OUT )

  set(COMM_ARRAY ${ARGN})

  if (COMM_ARRAY)
    set(ADD_SERIAL_FEATURE OFF)
    set(ADD_MPI_FEATURE OFF)
    foreach(COMM ${COMM_ARRAY})
      if(${COMM} STREQUAL "serial")
        set(ADD_SERIAL_FEATURE ON)
      elseif (${COMM} STREQUAL "mpi")
        set(ADD_MPI_FEATURE ON)
      else()
        message(SEND_ERROR "Error, the COMM value '${COMM}' is not valid."
          "  Only 'mpi' and 'serial' are allowed.")
      endif()
    endforeach()
  else()
    set(ADD_MPI_FEATURE ON)
    set(ADD_SERIAL_FEATURE ON)
  endif()

  if (TPL_ENABLE_MPI)
    set(ADD_SERIAL_FEATURE OFF)
    if (NOT ADD_MPI_FEATURE)
      set(EXCLUDED_FEATURE TRUE)
    endif()
  else()
    set(ADD_MPI_FEATURE OFF)
    if (NOT ADD_SERIAL_FEATURE)
      set(EXCLUDED_FEATURE TRUE)
    endif()
  endif()

  if (TEST_NAME AND EXCLUDED_FEATURE)
    message_wrapper(
      "-- ${TEST_NAME}: NOT added test because TPL_ENABLE_MPI='${TPL_ENABLE_MPI}' and COMM='${COMM_ARRAY}'!"
      )
  endif()

  set(${ADD_SERIAL_FEATURE_OUT} ${ADD_SERIAL_FEATURE} PARENT_SCOPE)
  set(${ADD_MPI_FEATURE_OUT} ${ADD_MPI_FEATURE}  PARENT_SCOPE)

endfunction()


function(tribits_create_name_from_current_source_directory  directoryNameOut)
    set(directoryName "")

    #Get the unique part of the path for this test directory
    string(REGEX REPLACE ${PACKAGE_SOURCE_DIR} "" unique_dir_path
      ${CMAKE_CURRENT_SOURCE_DIR})

    #strip off the preceding "/"
    string(LENGTH ${unique_dir_path} udp_length)
    math(EXPR last_index "${udp_length}-1")
    string(SUBSTRING ${unique_dir_path} 1 ${last_index} unique_dir_path)

    # Make the name acceptable for filesystems. This may need to be made
    # compatible with windows since they use a "\" instead of a "/" for
    # directory delimiters. I'm not sure how this will react if we encounter a
    # directory name with a space in it.
    string(REGEX REPLACE "/" "_" directoryName "${unique_dir_path}")

    set(${directoryNameOut} "${directoryName}" PARENT_SCOPE)
endfunction()

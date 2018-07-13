# -*- mode: cmake -*-
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
# Based on version from MSTK project which has the following:
# ACKNOWLEDGEMENT:
# This capability was developed for the Amanzi code base under the DOE/ASCEM
# project and used here in the MSTK project owing to its open source BSD license
# Many thanks to the Amanzi build system developers 
#

# 
# Print variable
#
#
# ADD_IMPORTED_LIBRARY(name [SHARED | STATIC]
#                      LOCATION <path>
#                      [ LINK_LANGUAGES <lang1> <lang2> <lang3> ... ]
#                      [ LINK_INTERFACE_LIBRARIES <lib1> <lib2> ... ]
#                     )
#                    
#                      

include(CMakeParseArguments)
function(ADD_IMPORTED_LIBRARY target_name)

  set(_options SHARED STATIC)
  set(_oneValueArgs LOCATION)
  set(_multiValueArgs LINK_LANGUAGES LINK_INTERFACE_LIBRARIES)

  cmake_parse_arguments(PARSE "${_options}" "${_oneValueArgs}" "${_multiValueArgs}" ${ARGN} ) 

  # --- Check what has been passed in

  # SHARED and STATIC can not be set at the same time
  if ( "${PARSE_STATIC}" AND "${PARSE_SHARED}" )
    message(FATAL_ERROR "Can not specify imported library as shared and static.")
  endif()

  # Require a location
  if ( NOT PARSE_LOCATION )
    message(FATAL_ERROR "Must specify a location to define an imported library target.")
  endif()

  # Check to see if name already exists as a target
  if ( NOT TARGET "${target_name}" )

    # --- Set the library type
    set(lib_type UNKNOWN)
    if(PARSE_STATIC)
      set(lib_type STATIC)
    endif()
    if(PARSE_SHARED)
      set(lib_type SHARED)
    endif()

    # --- Add the library target 
    add_library(${target_name} ${lib_type} IMPORTED)

    # --- Update the global property that tracks imported targets
    set(prop_name IMPORTED_${lib_type}_LIBRARIES)
    get_property(prop_value GLOBAL PROPERTY ${prop_name})
    set_property(GLOBAL PROPERTY ${prop_name} ${prop_value} ${target_name})

    # --- Set the properties
    set_target_properties(${target_name} PROPERTIES
                          IMPORTED_LOCATION ${PARSE_LOCATION})
    if ( PARSE_LINK_LANGUAGES )
      set_target_properties(${target_name} PROPERTIES
                          IMPORTED_LINK_INTERFACE_LANGUAGES "${PARSE_LINK_LANGUAGES}")
    endif()
    if ( PARSE_LINK_INTERFACE_LIBRARIES )
      set_target_properties(${target_name} PROPERTIES
                            IMPORTED_LINK_INTERFACE_LIBRARIES "${PARSE_LINK_INTERFACE_LIBRARIES}")
    endif()
     
  endif()
endfunction(ADD_IMPORTED_LIBRARY)


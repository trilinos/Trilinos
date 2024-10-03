# -*- mode: cmake -*-
# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
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
# add_imported_library(name [SHARED | STATIC]
#                      LOCATION <path>
#                      [ LINK_LANGUAGES <lang1> <lang2> <lang3> ... ]
#                      [ LINK_INTERFACE_LIBRARIES <lib1> <lib2> ... ]
#                     )
#                    
#                      

include(CMakeParseArguments)
function(add_imported_library target_name)

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
endfunction()


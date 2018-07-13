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
# Based on version from MSTK which is from Amanzi open source code -
# https://software.lanl.gov/ascem/trac)
#
#       ADD_PACKAGE_DEPENDENCY(<PACKNAME> DEPENDS_ON <req_pack>)
#

# CMake module
include(CMakeParseArguments)

# MSTK modules
include(ParseLibraryList)

function(ADD_PACKAGE_DEPENDENCY)

    # Macro: _print_usage
    macro(_print_usage)
        message("\nADD_PACKAGE_DEPENDENCY(<target_package> DEPENDS_ON <req_package>)\n"
                " Add req_package to target_package dependencies.\n")
    endmacro(_print_usage)

    # Parse the arguments
    set(_options "")
    set(_oneValue "DEPENDS_ON")
    set(_multiValue "")
    cmake_parse_arguments(ADD_PACK "${_options}" "${_oneValue}" "${_multiValue}" ${ARGN})
     
    # Define the target package name
    list(GET ADD_PACK_UNPARSED_ARGUMENTS 0 target_package)
    if ( NOT target_package )
        _print_usage()
        message(FATAL_ERROR "Must define a target_package")
    endif()

    # Define the required package
    set(req_package "")
    if(ADD_PACK_DEPENDS_ON)
        set(req_package ${ADD_PACK_DEPENDS_ON})
    else()
        _print_usage()
        message(FATAL_ERROR "Must define a required package")
    endif()    

    # Find required package
    message(STATUS "${target_package} depends on ${req_package}")
    message(STATUS "Updating ${target_package}_LIBRARIES and ${target_package}_INCLUDE_DIRS")
    find_package(${req_package} REQUIRED)
    if( ${req_package}_LIBRARIES AND ${req_package}_INCLUDE_DIRS )

        # Add the include paths 
        set(_save_inc_list ${${target_package}_INCLUDE_DIRS})
        list(APPEND _save_inc_list ${${req_package}_INCLUDE_DIRS})
        set(${target_package}_INCLUDE_DIRS ${_save_inc_list} PARENT_SCOPE)
        list(REMOVE_DUPLICATES ${target_package}_INCLUDE_DIRS)

        # Add the libraries....this can be tricky because some packages
        # *cough* HDF5 *cough* return a list with keywords debug, optimized
        # general in the list. These keywords are flags that are used when
        # CMAKE_BUILD_TYPE is set. Need to construct the LIBRARIES carefully
        # if these are present.
        parse_library_list(${${target_package}_LIBRARIES}
                           FOUND     target_libs_split
                           DEBUG     target_debug_libs
                           OPT       target_opt_libs
                           GENERAL   target_gen_libs)

        parse_library_list(${${req_package}_LIBRARIES}
                           FOUND     req_libs_split
                           DEBUG     req_debug_libs
                           OPT       req_opt_libs
                           GENERAL   req_gen_libs)
        
        # _save_lib_list tmp storage
        set(_save_lib_list "")
        if ( ${target_libs_split} OR ${req_libs_split} )

            # Define the parsed lists if the orginal list did not contain keywords
            if ( NOT ${target_libs_split} )
                set(target_debug_libs ${${target_package}_LIBRARIES})
                set(target_opt_libs ${${target_package}_LIBRARIES})
                set(target_general_libs ${${target_package}_LIBRARIES})
            endif()    

            if ( NOT ${req_libs_split} )
                set(req_debug_libs ${${req_package}_LIBRARIES})
                set(req_opt_libs ${${req_package}_LIBRARIES})
                set(req_general_libs ${${req_package}_LIBRARIES})
            endif()    

            # Define each type and store in tmp lists removing duplicates
            set(_save_debug_list "")
            set(_save_opt_list "")
            set(_save_gen_list "")
            if ( target_debug_libs OR req_debug_libs )
                set(_save_debug_list "${target_debug_libs}" "${req_debug_libs}")
                list(REMOVE_DUPLICATES _save_debug_list)
            endif()    
            if ( target_opt_libs OR req_opt_libs )
                set(_save_opt_list "${target_opt_libs}" "${req_opt_libs}")
                list(REMOVE_DUPLICATES _save_opt_list)
            endif()    
            if ( target_gen_libs OR req_gen_libs )
                set(_save_gen_list "${target_gen_libs}" "${req_gen_libs}")
                list(REMOVE_DUPLICATES _save_gen_list)
            endif()    

            # Now build the _save_lib_list with the keywords
            if(_save_debug_list)
                list(APPEND _save_lib_list "debug")
                list(APPEND _save_lib_list "${_save_debug_list}")
            endif()    
            
            if(_save_opt_list)
                list(APPEND _save_lib_list "optimized")
                list(APPEND _save_lib_list "${_save_opt_list}")
            endif()    

            if(_save_gen_list)
                list(APPEND _save_lib_list "general")
                list(APPEND _save_lib_list "${_save_gen_list}")
            endif()    

        else()

            #  Neither list has keywords
            set(_save_lib_list "${${target_package}_LIBRARIES}" "${${req_package}_LIBRARIES}")
            list(REMOVE_DUPLICATES _save_lib_list)

        endif()    

       
        set(${target_package}_LIBRARIES ${_save_lib_list} PARENT_SCOPE)

    endif()    

endfunction(ADD_PACKAGE_DEPENDENCY)

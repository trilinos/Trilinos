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
# Based on the MSTK CGNS Find Module
#
# Usage:
#    Control the search through CGNS_DIR or setting environment variable
#    CGNS_ROOT to the CGNS installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    CGNS_FOUND            (BOOL)       Flag indicating if CGNS was found
#    CGNS_INCLUDE_DIR      (PATH)       Path to the CGNS include file
#    CGNS_INCLUDE_DIRS     (LIST)       List of all required include files
#    CGNS_LIBRARY_DIR      (PATH)       Path to the CGNS library
#    CGNS_LIBRARY          (FILE)       CGNS library
#    CGNS_LIBRARIES        (LIST)       List of all required CGNS libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# MSTK CMake functions see <root>/cmake/modules for source
include(AddPackageDependency)

# If CGNS_ROOT was defined in the environment, use it.
# Definition from the command line will take precedence.
if (NOT CGNS_ROOT AND NOT $ENV{CGNS_ROOT} STREQUAL "")
  set(CGNS_ROOT $ENV{CGNS_ROOT})
endif()

# CGNS_DIR is DEPRECATED WARN THE USER if it is set
if (NOT CGNS_ROOT AND CGNS_DIR )
  message(WARNING "The configuration parameter CGNS_DIR is deprecated."
                  " Please use CGNS_ROOT instead to define the CGNS installation")
  set(CGNS_ROOT ${CGNS_DIR})
endif()  

# Add the usual paths for searching using the CGNS_ROOT variable
if (CGNS_ROOT)
  list(APPEND _cgns_INCLUDE_SEARCH_DIRS 
              ${CGNS_ROOT}/include
              ${CGNS_ROOT})
 
  list(APPEND _cgns_LIBRARY_SEARCH_DIRS 
              ${CGNS_ROOT}/lib
              ${CGNS_ROOT})
  
            list(APPEND _cgns_BINARY_SEARCH_DIRS 
              ${CGNS_ROOT}/bin
              ${CGNS_ROOT})
endif()
 
if ( CGNS_LIBRARIES AND CGNS_INCLUDE_DIRS )

    # Do nothing. Variables are set. No need to search again

else(CGNS_LIBRARIES AND CGNS_INCLUDE_DIRS)

    # Cache variables
    if(CGNS_ROOT)
        set(CGNS_ROOT "${CGNS_ROOT}" CACHE PATH "Path to search for CGNS include and library files")
    endif()

    if(CGNS_INCLUDE_DIR)
        set(CGNS_INCLUDE_DIR "${CGNS_INCLUDE_DIR}" CACHE PATH "Path to search for CGNS include files")
    endif()

    if(CGNS_LIBRARY_DIR)
        set(CGNS_LIBRARY_DIR "${CGNS_LIBRARY_DIR}" CACHE PATH "Path to search for CGNS library files")
    endif()

    
    # Search for include files
    # Search order preference:
    #  (1) CGNS_INCLUDE_DIR - check existence of path AND if the include files exist
    #  (2) CGNS_ROOT/<include>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(cgns_inc_names "cgnslib.h")
    if (CGNS_INCLUDE_DIR)

        if (EXISTS "${CGNS_INCLUDE_DIR}")

            find_path(cdf_test_include_path
                      NAMES ${cgns_inc_names}
                      HINTS ${CGNS_INCLUDE_DIR}
                      NO_DEFAULT_PATH)
            if(NOT cdf_test_include_path)
                message(SEND_ERROR "Can not locate ${cgns_inc_names} in ${CGNS_INCLUDE_DIR}")
            endif()
            set(CGNS_INCLUDE_DIR "${cdf_test_include_path}")

        else()
            message(SEND_ERROR "CGNS_INCLUDE_DIR=${CGNS_INCLUDE_DIR} does not exist")
            set(CGNS_INCLUDE_DIR "CGNS_INCLUDE_DIR-NOTFOUND")
        endif()

    else() 

        set(cgns_inc_suffixes "include")
        if(CGNS_ROOT)

            if (EXISTS "${CGNS_ROOT}" )

                find_path(CGNS_INCLUDE_DIR
                          NAMES ${cgns_inc_names}
                          HINTS ${CGNS_ROOT}/include
                          PATH_SUFFIXES ${cgns_inc_suffixes}
                          NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "CGNS_ROOT=${CGNS_ROOT} does not exist")
                 set(CGNS_INCLUDE_DIR "CGNS_INCLUDE_DIR-NOTFOUND")
            endif()    


        else()

            find_path(CGNS_INCLUDE_DIR
                      NAMES ${cgns_inc_names}
                      PATH_SUFFIXES ${cgns_inc_suffixes})

        endif()

    endif()


    if ( NOT CGNS_INCLUDE_DIR )
        message(SEND_ERROR "Can not locate CGNS include directory")
    endif()

    # Search for libraries 
    # Search order preference:
    #  (1) CGNS_LIBRARY_DIR - check existence of path AND if the include files exist
    #  (2) CGNS_ROOT/<lib,Lib>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    if (CGNS_LIBRARY_DIR)

        if (EXISTS "${CGNS_LIBRARY_DIR}")

            find_library(CGNS_LIBRARY
                         NAMES cgns
                         HINTS ${CGNS_LIBRARY_DIR}
                         NO_DEFAULT_PATH)

        else()
            message(SEND_ERROR "CGNS_LIBRARY_DIR=${CGNS_LIBRARY_DIR} does not exist")
            set(CGNS_LIBRARY "CGNS_LIBRARY-NOTFOUND")
        endif()

    else() 

        if(CGNS_ROOT)

            if (EXISTS "${CGNS_ROOT}" )

                find_library(CGNS_LIBRARY
                             NAMES cgns
                             HINTS ${CGNS_ROOT}
                             PATH_SUFFIXES "lib" "Lib"
                             NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "CGNS_ROOT=${CGNS_ROOT} does not exist")
                 set(CGNS_LIBRARY "CGNS_LIBRARY-NOTFOUND")
            endif()    


        else()

            find_library(CGNS_LIBRARY
                         NAMES cgns
                         PATH_SUFFIXES ${cgns_lib_suffixes})
            
        endif()

    endif()

    if ( NOT CGNS_LIBRARY )
        message(SEND_ERROR "Can not locate CGNS library")
    endif()    
    
    # Define the LIBRARIES and INCLUDE_DORS
    set(CGNS_INCLUDE_DIRS ${CGNS_INCLUDE_DIR})
    set(CGNS_LIBRARIES    ${CGNS_CXX_LIBRARY} ${CGNS_LIBRARY})

    # Need to find the CGNS config script to check for HDF5
    set(cgns_config_h "${CGNS_INCLUDE_DIR}/cgnsconfig.h" )
    if (EXISTS ${cgns_config_h})
      file(STRINGS "${cgns_config_h}" cg_build_hdf5_string REGEX "^#define CG_BUILD_HDF5")
      string(REGEX REPLACE "[^0-9]" "" cg_build_hdf5 "${cg_build_hdf5_string}")
      if ( cg_build_hdf5 EQUAL 51 ) # Kluge: define is 1, but the 5 comes from hdf5
          message(STATUS "CGNS requires HDF5")
          add_package_dependency(CGNS DEPENDS_ON HDF5)
      endif()
    else()
      message(STATUS "CGNS does not have cgnsconfig.h; assuming CGNS depends on HDF5")     
      add_package_dependency(CGNS DEPENDS_ON HDF5)
    endif()

    # Need to find the CGNS types include to check for SCOPING
    set(cgns_types_h "${CGNS_INCLUDE_DIR}/cgnstypes.h" )
    if (EXISTS ${cgns_types_h})
      file(STRINGS "${cgns_types_h}" cg_scoping_string REGEX "^#define CG_BUILD_SCOPE")
      string(REGEX REPLACE "[^0-9]" "" cg_scoping "${cg_scoping_string}")
      if ( cg_scoping EQUAL 1 )
          message(STATUS "CGNS Scoping Enabled as Required")
      else()
          message(SEND_ERROR "CGNS Scoping *Not* Enabled as Required. Rebuild CGNS library with CGNS_ENABLE_SCOPING defined.")
      endif()
    else()
      message(SEND_ERROR "CGNS: Could not find cgnstypes.h")
    endif()
endif(CGNS_LIBRARIES AND CGNS_INCLUDE_DIRS )    

# Search for CGNS tools
set(_cgns_TOOLS cgnscheck cgnsdiff cgnslist cgnscompress cgnsconvert cgnsnames cgnsupdate)
set(CGNS_TOOLS_FOUND)
foreach( tool ${_cgns_TOOLS})
  string(TOUPPER "${tool}" tool_uc)
  set(_cgns_VAR_NAME CGNS_${tool_uc}_BINARY)
  find_program(${_cgns_VAR_NAME}
               ${tool}
               HINTS ${_cgns_BINARY_SEARCH_DIRS}
               ${_cgns_FIND_OPTIONS})
  if (${_cgns_VAR_NAME})
    list(APPEND CGNS_TOOLS_FOUND ${tool})
  endif()
endforeach()

# Send useful message if everything is found
find_package_handle_standard_args(CGNS DEFAULT_MSG
                                        CGNS_LIBRARIES
                                        CGNS_INCLUDE_DIRS)

# find_package)handle)standard_args should set CGNS_FOUND but it does not!
if ( CGNS_LIBRARIES AND CGNS_INCLUDE_DIRS)
    set(CGNS_FOUND TRUE)
else()
    set(CGNS_FOUND FALSE)
endif()

# --- Provide a summary of what the module found
if ( NOT CGNS_FIND_QUIETLY )

  # Create a not found list
  message(STATUS "\tCGNS_INCLUDE_DIRS      = ${CGNS_INCLUDE_DIRS}")
  message(STATUS "\tCGNS_LIBRARIES         = ${CGNS_LIBRARIES}")
  message(STATUS "\tCGNS_TOOLS_FOUND       = ${CGNS_TOOLS_FOUND}")

endif()
# For compatibility with TriBITS:
set(DOCSTR "List of semi-colon separated paths to look for the TPL CGNS")

set(TPL_CGNS_LIBRARIES ${CGNS_LIBRARIES} CACHE PATH ${DOCSTR})
set(TPL_CGNS_INCLUDE_DIRS ${CGNS_INCLUDE_DIRS} CACHE PATH ${DOCSTR})

mark_as_advanced(
  CGNS_INCLUDE_DIR
  CGNS_INCLUDE_DIRS
  CGNS_LIBRARY
  CGNS_CXX_LIBRARY
  CGNS_LIBRARIES
  CGNS_LIBRARY_DIR
)

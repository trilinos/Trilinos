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
# Based on the MSTK PNetCDF Find Module
#
# Usage:
#    Control the search through PNetCDF_DIR or setting environment variable
#    PNetCDF_ROOT to the PNetCDF installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    PNetCDF_FOUND            (BOOL)       Flag indicating if PNetCDF was found
#    PNetCDF_INCLUDE_DIR      (PATH)       Path to the PNetCDF include file
#    PNetCDF_INCLUDE_DIRS     (LIST)       List of all required include files
#    PNetCDF_LIBRARY_DIR      (PATH)       Path to the PNetCDF library
#    PNetCDF_LIBRARY          (FILE)       PNetCDF library
#    PNetCDF_LIBRARIES        (LIST)       List of all required PNetCDF libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# MSTK CMake functions see <root>/cmake/modules for source
include(AddPackageDependency)

if ( PNetCDF_LIBRARIES AND PNetCDF_INCLUDE_DIRS )

    # Do nothing. Variables are set. No need to search again

else(PNetCDF_LIBRARIES AND PNetCDF_INCLUDE_DIRS)

    # If NetCDF_ROOT was defined in the environment, use it.
    # Definition from the command line will take precedence.
    if (NOT PNetCDF_ROOT AND NOT $ENV{PNetCDF_ROOT} STREQUAL "")
      set(PNetCDF_ROOT $ENV{PNetCDF_ROOT})
    endif()

    # PNetCDF_DIR is DEPRECATED WARN THE USER if it is set
    if (NOT PNetCDF_ROOT AND PNetCDF_DIR )
      message(WARNING "The configuration parameter PNetCDF_DIR is deprecated."
                      " Please use PNetCDF_ROOT instead to define the NetCDF installation")
      set(PNetCDF_ROOT ${PNetCDF_DIR})
    endif()  

    # Cache variables
    if(PNetCDF_ROOT)
        set(PNetCDF_ROOT "${PNetCDF_ROOT}" CACHE PATH "Path to search for PNetCDF include and library files")
    endif()

    if(PNetCDF_INCLUDE_DIR)
        set(PNetCDF_INCLUDE_DIR "${PNetCDF_INCLUDE_DIR}" CACHE PATH "Path to search for PNetCDF include files")
    endif()

    if(PNetCDF_LIBRARY_DIR)
        set(PNetCDF_LIBRARY_DIR "${PNetCDF_LIBRARY_DIR}" CACHE PATH "Path to search for PNetCDF library files")
    endif()

    
    # Search for include files
    # Search order preference:
    #  (1) PNetCDF_INCLUDE_DIR - check existence of path AND if the include files exist
    #  (2) PNetCDF_ROOT/<include>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(pnetcdf_inc_names "pnetcdf.h")
    if (PNetCDF_INCLUDE_DIR)

        if (EXISTS "${PNetCDF_INCLUDE_DIR}")

            find_path(cdf_test_include_path
                      NAMES ${pnetcdf_inc_names}
                      HINTS ${PNetCDF_INCLUDE_DIR}
                      NO_DEFAULT_PATH)
            if(NOT cdf_test_include_path)
                message(SEND_ERROR "Can not locate ${pnetcdf_inc_names} in ${PNetCDF_INCLUDE_DIR}")
            endif()
            set(PNetCDF_INCLUDE_DIR "${cdf_test_include_path}")

        else()
            message(SEND_ERROR "PNetCDF_INCLUDE_DIR=${PNetCDF_INCLUDE_DIR} does not exist")
            set(PNetCDF_INCLUDE_DIR "PNetCDF_INCLUDE_DIR-NOTFOUND")
        endif()

    else() 

        set(pnetcdf_inc_suffixes "include")
        if(PNetCDF_ROOT)

            if (EXISTS "${PNetCDF_ROOT}" )

                find_path(PNetCDF_INCLUDE_DIR
                          NAMES ${pnetcdf_inc_names}
                          HINTS ${PNetCDF_ROOT}/include
                          PATH_SUFFIXES ${pnetcdf_inc_suffixes}
                          NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "PNetCDF_ROOT=${PNetCDF_ROOT} does not exist")
                 set(PNetCDF_INCLUDE_DIR "PNetCDF_INCLUDE_DIR-NOTFOUND")
            endif()    


        else()

            find_path(PNetCDF_INCLUDE_DIR
                      NAMES ${pnetcdf_inc_names}
                      PATH_SUFFIXES ${pnetcdf_inc_suffixes})

        endif()

    endif()


    if ( NOT PNetCDF_INCLUDE_DIR )
        message(SEND_ERROR "Can not locate PNetCDF include directory")
    endif()

    # Large dimension and parallel check here
    if ( PNetCDF_INCLUDE_DIR ) 
       
        set(pnetcdf_h "${PNetCDF_INCLUDE_DIR}/pnetcdf.h" )
        message(STATUS "PNetCDF include file ${pnetcdf_h} will be searched for define values")

        file(STRINGS "${pnetcdf_h}" pnetcdf_max_dims_string REGEX "^#define NC_MAX_DIMS")
        string(REGEX REPLACE "[^0-9]" "" pnetcdf_max_dims "${pnetcdf_max_dims_string}")

        file(STRINGS "${pnetcdf_h}" pnetcdf_max_vars_string REGEX "^#define NC_MAX_VARS")
        string(REGEX REPLACE "[^0-9]" "" pnetcdf_max_vars "${pnetcdf_max_vars_string}")

        if ( 
             ( (pnetcdf_max_dims EQUAL 65536)  OR (pnetcdf_max_dims GREATER 65536) ) AND
             ( (pnetcdf_max_vars EQUAL 524288) OR (pnetcdf_max_vars GREATER 524288) )
            )
            set(PNetCDF_LARGE_DIMS TRUE)
        else()
            message(WARNING "WARNING: The PNetCDF found in ${PNetCDF_ROOT} does not have the correct NC_MAX_DIMS and NC_MAX_VARS. "
                             "It may not be compatible with Exodus. See NetCDF-Mapping.md for details\n" )
            set(PNetCDF_LARGE_DIMS FALSE)
        endif()

    endif()    

    # Search for libraries 
    # Search order preference:
    #  (1) PNetCDF_LIBRARY_DIR - check existence of path AND if the include files exist
    #  (2) PNetCDF_ROOT/<lib,Lib>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information
    #
    if (PNetCDF_LIBRARY_DIR)

        if (EXISTS "${PNetCDF_LIBRARY_DIR}")

            find_library(PNetCDF_LIBRARY
                         NAMES pnetcdf
                         HINTS ${PNetCDF_LIBRARY_DIR}
                         NO_DEFAULT_PATH)

        else()
            message(SEND_ERROR "PNetCDF_LIBRARY_DIR=${PNetCDF_LIBRARY_DIR} does not exist")
            set(PNetCDF_LIBRARY "PNetCDF_LIBRARY-NOTFOUND")
        endif()

    else() 

        if(PNetCDF_ROOT)

            if (EXISTS "${PNetCDF_ROOT}" )

                find_library(PNetCDF_LIBRARY
                             NAMES pnetcdf
                             HINTS ${PNetCDF_ROOT}
                             PATH_SUFFIXES "lib" "Lib"
                             NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "PNetCDF_ROOT=${PNetCDF_ROOT} does not exist")
                 set(PNetCDF_LIBRARY "PNetCDF_LIBRARY-NOTFOUND")
            endif()    


        else()

            find_library(PNetCDF_LIBRARY
                         NAMES pnetcdf
                         PATH_SUFFIXES ${pnetcdf_lib_suffixes})
            
        endif()

    endif()

    if ( NOT PNetCDF_LIBRARY )
        message(SEND_ERROR "Can not locate PNetCDF library")
    endif()    
    
    # Define the LIBRARIES and INCLUDE_DORS
    set(PNetCDF_INCLUDE_DIRS ${PNetCDF_INCLUDE_DIR})
    set(PNetCDF_LIBRARIES    ${PNetCDF_CXX_LIBRARY} ${PNetCDF_LIBRARY})

endif(PNetCDF_LIBRARIES AND PNetCDF_INCLUDE_DIRS )    

# Send useful message if everything is found
find_package_handle_standard_args(PNetCDF DEFAULT_MSG
                                        PNetCDF_LIBRARIES
                                        PNetCDF_INCLUDE_DIRS)

# find_package)handle)standard_args should set PNetCDF_FOUND but it does not!
if ( PNetCDF_LIBRARIES AND PNetCDF_INCLUDE_DIRS)
    set(PNetCDF_FOUND TRUE)
else()
    set(PNetCDF_FOUND FALSE)
endif()

# --- Provide a summary of what the module found
if ( NOT PNetCDF_FIND_QUIETLY )

  # Create a not found list
  message(STATUS "\tPNetCDF_INCLUDE_DIRS      = ${PNetCDF_INCLUDE_DIRS}")
  message(STATUS "\tPNetCDF_LIBRARIES         = ${PNetCDF_LIBRARIES}")

endif()
# For compatibility with TriBITS:
set(DOCSTR "List of semi-colon separated paths to look for the TPL PNetCDF")

set(TPL_PNetCDF_LIBRARIES ${PNetCDF_LIBRARIES} CACHE PATH ${DOCSTR})
set(TPL_PNetCDF_INCLUDE_DIRS ${PNetCDF_INCLUDE_DIRS} CACHE PATH ${DOCSTR})

mark_as_advanced(
  PNetCDF_INCLUDE_DIR
  PNetCDF_INCLUDE_DIRS
  PNetCDF_LIBRARY
  PNetCDF_CXX_LIBRARY
  PNetCDF_LIBRARIES
  PNetCDF_LIBRARY_DIR
)

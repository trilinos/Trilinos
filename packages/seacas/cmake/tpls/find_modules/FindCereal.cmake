# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
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
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER
#
# Based on the CGNS Find Module in this repository
#
# Usage:
#    Control the search through Cereal_DIR or setting environment variable
#    Cereal_ROOT to the Cereal installation prefix.
#
#    This module does not search default paths!
#
#    Following variables are set:
#    Cereal_FOUND            (BOOL)       Flag indicating if Cereal was found
#    Cereal_INCLUDE_DIR      (PATH)       Path to the Cereal include file
#    Cereal_INCLUDE_DIRS     (LIST)       List of all required include files
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# MSTK CMake functions see <root>/cmake/modules for source
include(AddPackageDependency)

# If Cereal_ROOT was defined in the environment, use it.
# Definition from the command line will take precedence.
if (NOT Cereal_ROOT AND NOT $ENV{Cereal_ROOT} STREQUAL "")
  set(Cereal_ROOT $ENV{Cereal_ROOT})
endif()

# Cereal_DIR is DEPRECATED WARN THE USER if it is set
if (NOT Cereal_ROOT AND Cereal_DIR )
  message(WARNING "The configuration parameter Cereal_DIR is deprecated."
    " Please use Cereal_ROOT instead to define the Cereal installation")
  set(Cereal_ROOT ${Cereal_DIR})
endif()

# Add the usual paths for searching using the Cereal_ROOT variable
if (Cereal_ROOT)
  list(APPEND _cereal_INCLUDE_SEARCH_DIRS
    ${Cereal_ROOT}/include
    ${Cereal_ROOT})
endif()

if ( Cereal_INCLUDE_DIRS )

  # Do nothing. Variables are set. No need to search again

else( Cereal_INCLUDE_DIRS )

  # Cache variables
  if(Cereal_ROOT)
    set(Cereal_ROOT "${Cereal_ROOT}" CACHE PATH "Path to search for Cereal include and library files")
  endif()

  if(Cereal_INCLUDE_DIR)
    set(Cereal_INCLUDE_DIR "${Cereal_INCLUDE_DIR}" CACHE PATH "Path to search for Cereal include files")
  endif()

  # Search for include files
  # Search order preference:
  #  (1) Cereal_INCLUDE_DIR - check existence of path AND if the include files exist
  #  (2) Cereal_ROOT/<include>
  #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
  #
  set(cereal_inc_names "cereal/cereal.hpp")
  if (Cereal_INCLUDE_DIR)

    if (EXISTS "${Cereal_INCLUDE_DIR}")

      find_path(cdf_test_include_path
        NAMES ${cereal_inc_names}
        HINTS ${Cereal_ROOT}/include
        NO_DEFAULT_PATH)
      if(NOT cdf_test_include_path)
        message(SEND_ERROR "Can not locate ${cereal_inc_names} in ${Cereal_INCLUDE_DIR}")
      endif()
      set(Cereal_INCLUDE_DIR "${cdf_test_include_path}")

    else()
      message(SEND_ERROR "Cereal_INCLUDE_DIR=${Cereal_INCLUDE_DIR} does not exist")
      set(Cereal_INCLUDE_DIR "Cereal_INCLUDE_DIR-NOTFOUND")
    endif()

  else()

    set(cereal_inc_suffixes "include")
    if(Cereal_ROOT)

      if (EXISTS "${Cereal_ROOT}" )

        find_path(Cereal_INCLUDE_DIR
          NAMES ${cereal_inc_names}
          HINTS ${Cereal_ROOT}/include
          PATH_SUFFIXES ${cereal_inc_suffixes}
          NO_DEFAULT_PATH)

      else()
        message(SEND_ERROR "Cereal_ROOT=${Cereal_ROOT} does not exist")
        set(Cereal_INCLUDE_DIR "Cereal_INCLUDE_DIR-NOTFOUND")
      endif()


    else()

      find_path(Cereal_INCLUDE_DIR
        NAMES ${cereal_inc_names}
        PATH_SUFFIXES ${cereal_inc_suffixes})

    endif()

  endif()


  if ( NOT Cereal_INCLUDE_DIR )
    message(SEND_ERROR "Can not locate Cereal include directory")
  endif()

  # Define the INCLUDE_DIRS
  set(Cereal_INCLUDE_DIRS ${Cereal_INCLUDE_DIR})

endif( Cereal_INCLUDE_DIRS )


# Send useful message if everything is found
find_package_handle_standard_args(Cereal DEFAULT_MSG
  Cereal_INCLUDE_DIRS)

# find_package)handle)standard_args should set Cereal_FOUND but it does not!
if ( Cereal_INCLUDE_DIRS)
  set(Cereal_FOUND TRUE)
else()
  set(Cereal_FOUND FALSE)
endif()

# --- Provide a summary of what the module found
if ( NOT Cereal_FIND_QUIETLY )

  # Create a not found list
  message(STATUS "\tCereal_FOUND             = ${Cereal_FOUND}")
  message(STATUS "\tCereal_INCLUDE_DIR       = ${Cereal_INCLUDE_DIR}")
  message(STATUS "\tCereal_INCLUDE_DIRS      = ${Cereal_INCLUDE_DIRS}")

endif()
# For compatibility with TriBITS:
SET(DOCSTR "List of semi-colon separated paths to look for the TPL Cereal")

set(TPL_Cereal_INCLUDE_DIRS ${Cereal_INCLUDE_DIRS} CACHE PATH ${DOCSTR})

mark_as_advanced(
  Cereal_INCLUDE_DIR
  Cereal_INCLUDE_DIRS
  )

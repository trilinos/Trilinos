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

INCLUDE(AssertDefined)
INCLUDE(AppendSet)
INCLUDE(SetNotFound)

# Search for all MUMPS libraries if they are not explicitly stated.
IF (NOT MUMPS_LIBRARY_NAMES OR NOT TPL_MUMPS_LIBRARIES)

  IF (MUMPS_LIBRARY_DIRS)
    SET(PATHS_ARG PATHS ${MUMPS_LIBRARY_DIRS})
  ELSE()
    SET(PATHS_ARG PATHS)
  ENDIF()

  SET(MUMPS_LIBRARIES
    smumps dmumps cmumps zmumps)
  SET(MUMPS_LIBRARIES_FOUND)
  SET(MUMPS_HEADERS)

  FOREACH(LIBNAME ${MUMPS_LIBRARIES})

    SET_NOTFOUND( _MUMPS_${LIBNAME}_LIBRARY )
    FIND_LIBRARY( _MUMPS_${LIBNAME}_LIBRARY
      NAMES ${LIBNAME}
      ${PATHS_ARG} NO_DEFAULT_PATH )
    FIND_LIBRARY( _MUMPS_${LIBNAME}_LIBRARY
      NAMES ${LIBNAME} )

    IF (_MUMPS_${LIBNAME}_LIBRARY)
      APPEND_SET(MUMPS_LIBRARIES_FOUND ${LIBNAME})
      APPEND_SET(MUMPS_HEADERS ${LIBNAME}_c.h)
    ENDIF()
 
  ENDFOREACH()

  #Add ordering and common library
  SET(MUMPS_COMMON_LIBRARIES
    pord mumps_common)
  FOREACH(LIBNAME ${MUMPS_COMMON_LIBRARIES})

    SET_NOTFOUND( _MUMPS_${LIBNAME}_LIBRARY )
    FIND_LIBRARY( _MUMPS_${LIBNAME}_LIBRARY
      NAMES ${LIBNAME}
      ${PATHS_ARG} NO_DEFAULT_PATH )
    FIND_LIBRARY( _MUMPS_${LIBNAME}_LIBRARY
      NAMES ${LIBNAME} )

    IF (_MUMPS_${LIBNAME}_LIBRARY)
      APPEND_SET(MUMPS_LIBRARIES_FOUND ${LIBNAME})
    ENDIF()
 
  ENDFOREACH()

  IF (NOT MUMPS_LIBRARIES_FOUND)
    MESSAGE(
      "-- ERROR: Did not find a lib in the lib set \"${MUMPS_LIBRARIES}\""
      " for the TPL 'MUMPS'!")
  ENDIF()

  TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( MUMPS
    REQUIRED_HEADERS ${MUMPS_HEADERS}
    REQUIRED_LIBS_NAMES ${MUMPS_LIBRARIES_FOUND}
  )

ELSE()

  # Check the given library names in the standard way.
  TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( MUMPS
    REQUIRED_HEADERS dmumps_c.h
    REQUIRED_LIBS_NAMES dmumps pord
  )

ENDIF()


# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
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

# This file gets included in the base-level CMakeLists.txt file to define
# Fortran name mangling.

IF (${PROJECT_NAME}_ENABLE_Fortran)
  INCLUDE(FortranMangling)
  FORTRAN_MANGLING()

  # Verify the selected combination of Fortran and C++ compilers.
  IF("${CMAKE_VERSION}" VERSION_GREATER 2.7.20090824 AND NOT ${PROJECT_NAME}_SKIP_FORTRANCINTERFACE_VERIFY_TEST)
    INCLUDE(FortranCInterface)
    FortranCInterface_VERIFY(CXX)
  ENDIF()
ENDIF()

IF (FC_FUNC_DEFAULT)

  SET(F77_FUNC_DEFAULT ${FC_FUNC_DEFAULT})
  SET(F77_FUNC__DEFAULT ${FC_FUNC__DEFAULT})
  # 2008/10/26: rabartl: ToDo: Above, we need to write
  # a different function to find out the right BLAS
  # name mangling automatically.  Given what the above
  # FORTRAN_MANGLING() function does, this should not
  # be too hard.

ELSE()
 
  IF(CYGWIN)
    SET(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
    SET(F77_FUNC__DEFAULT "(name,NAME) name ## __" )
  ELSEIF(WIN32)
    SET(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
    SET(F77_FUNC__DEFAULT "(name,NAME) NAME")
  ELSEIF(UNIX AND NOT APPLE)
    SET(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
    #SET(F77_FUNC__DEFAULT "(name,NAME) name ## __" )
    SET(F77_FUNC__DEFAULT "(name,NAME) name ## _" )
  ELSEIF(APPLE)
    SET(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
    SET(F77_FUNC__DEFAULT "(name,NAME) name ## __" )
  ELSE()
    MESSAGE(FATAL_ERROR "Error, could not determine fortran name mangling!")
  ENDIF()

ENDIF()

# Set options so that users can change these!

SET(F77_FUNC ${F77_FUNC_DEFAULT} CACHE STRING
  "Name mangling used to call Fortran 77 functions with no underscores in the name")
SET(F77_FUNC_ ${F77_FUNC__DEFAULT} CACHE STRING
  "Name mangling used to call Fortran 77 functions with at least one underscore in the name")

MARK_AS_ADVANCED(F77_FUNC)
MARK_AS_ADVANCED(F77_FUNC_)

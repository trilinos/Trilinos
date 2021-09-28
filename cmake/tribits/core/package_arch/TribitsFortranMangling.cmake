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

# This file gets included in the base-level CMakeLists.txt file to define
# Fortran name mangling.

if (${PROJECT_NAME}_ENABLE_CXX AND ${PROJECT_NAME}_ENABLE_Fortran)
  include(FortranMangling)
  fortran_mangling()

  # Verify the selected combination of Fortran and C++ compilers.
  if(NOT ${PROJECT_NAME}_SKIP_FORTRANCINTERFACE_VERIFY_TEST)
    include(FortranCInterface)
    fortrancinterface_verify(CXX)
  endif()
endif()

if (FC_FUNC_DEFAULT)

  set(F77_FUNC_DEFAULT ${FC_FUNC_DEFAULT})
  set(F77_FUNC__DEFAULT ${FC_FUNC__DEFAULT})
  # 2008/10/26: rabartl: ToDo: Above, we need to write
  # a different function to find out the right BLAS
  # name mangling automatically.  Given what the above
  # fortran_mangling() function does, this should not
  # be too hard.

else()

  if(CYGWIN)
    set(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
    set(F77_FUNC__DEFAULT "(name,NAME) name ## __" )
  elseif(WIN32)
    set(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
    set(F77_FUNC__DEFAULT "(name,NAME) NAME")
  elseif(UNIX AND NOT APPLE)
    set(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
    #set(F77_FUNC__DEFAULT "(name,NAME) name ## __" )
    set(F77_FUNC__DEFAULT "(name,NAME) name ## _" )
  elseif(APPLE)
    set(F77_FUNC_DEFAULT "(name,NAME) name ## _" )
    set(F77_FUNC__DEFAULT "(name,NAME) name ## __" )
  else()
    message(FATAL_ERROR "Error, could not determine fortran name mangling!")
  endif()

endif()

# Set options so that users can change these!

set(F77_FUNC ${F77_FUNC_DEFAULT} CACHE STRING
  "Name mangling used to call Fortran 77 functions with no underscores in the name")
set(F77_FUNC_ ${F77_FUNC__DEFAULT} CACHE STRING
  "Name mangling used to call Fortran 77 functions with at least one underscore in the name")

mark_as_advanced(F77_FUNC)
mark_as_advanced(F77_FUNC_)

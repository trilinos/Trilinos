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

# Function that defines variables describing the BLAS name mangling convention
#
# Sets the following outputs on success:
#
#  BLAS_FN_CASE
#    "UPPER" if names are translated to upper-case,
#    "LOWER" otherwise.
#
#  BLAS_FN_UNDERSCORE
#    "NO_UNDER" if nothing is appended, "UNDER" if
#    one underscore is appended, and "SECOND_UNDER"
#    if a function with an underscore has a second
#    appended.
#
# BLAS_FUNC_DEFAULT
#    The default mange mangling for Fortran functions
#    that do not contain an underscore.
#
#  The Fortran 2003 name binding facilities and ISO_C_BINDING module
#  should be preferred over cpp macro trickery whenever possible.
#
function(blas_mangling)

  if(NOT DEFINED BLAS_FN_CASE)

    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("BLAS_MANGLING: Testing name Mangling Schemes!\n")
    endif()

    find_file(_blascmakelists blasmangle/ ${CMAKE_MODULE_PATH})
    if (NOT _blascmakelists)
      message(STATUS "Error, the file blasmangle could not be found so we can not determine Fortran name mangling!")
      return()
    endif()

    set(_fcmangledir ${PROJECT_BINARY_DIR}/CMakeFiles/CMakeTmp/blasmangle)
    file(MAKE_DIRECTORY ${_fcmangledir})

    foreach(cdef LOWER UPPER)

      foreach(udef UNDER NO_UNDER)

        if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          message("BLAS_MANGLING: Testing ${cdef} ${udef}\n\n")
        endif()

        set(_fcmangledir_case "${_fcmangledir}/${cdef}/${udef}")
        file(MAKE_DIRECTORY "${_fcmangledir}/${cdef}")
        file(MAKE_DIRECTORY ${_fcmangledir_case})

        try_compile(_blasmngl ${_fcmangledir_case} ${_blascmakelists} blasmangle
          CMAKE_FLAGS "-DMANGLE_FLAGS:STRING=-DFC_FN_${cdef};-DFC_FN_${udef}"
          OUTPUT_VARIABLE _blasmngl_output
          )
        if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          message("${_blasmngl_output}\n\n")
        endif()

        if(_blasmngl)
          if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
            message("BLAS_MANGLING: Bingo!  ${cdef} ${udef} is the correct BLAS name mangling!\n")
          endif()
          set(BLAS_FN_CASE ${cdef} CACHE INTERNAL
            "Case used by Fortran functions" FORCE)
          set(BLAS_FN_UNDERSCORE ${udef} CACHE INTERNAL
            "Underscore convention used by Fortran functions" FORCE)
          break()
        endif()

      endforeach()

      if(_blasmngl)
        break()
      endif()

    endforeach()

    if(_blasmngl)
      message(STATUS "BLAS name mangling: ${BLAS_FN_CASE} ${BLAS_FN_UNDERSCORE}")
    else()
      message(STATUS "Warning, cannot automatically determine BLAS mangling.")
    endif()

  endif()

  if (BLAS_FN_CASE STREQUAL LOWER)
    set(BLAS_NAME_NAME name)
  elseif (BLAS_FN_CASE STREQUAL UPPER)
    set(BLAS_NAME_NAME NAME)
  endif()

  if (BLAS_FN_UNDERSCORE)
    if(BLAS_FN_UNDERSCORE STREQUAL "UNDER")
      set(BLAS_FUNC_DEFAULT "(name,NAME) ${BLAS_NAME_NAME} ## _" CACHE INTERNAL "")
    else()
      set(BLAS_FUNC_DEFAULT "(name,NAME) ${BLAS_NAME_NAME}" CACHE INTERNAL "")
    endif()
  endif()

endfunction()

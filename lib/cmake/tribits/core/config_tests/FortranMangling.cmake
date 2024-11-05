# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(GlobalSet)

# Function that defines variables describing the Fortran name mangling
# convention
#
# Sets the following outputs on success:
#
#  FC_FN_CASE
#    "UPPER" if names are translated to upper-case,
#    "LOWER" otherwise.
#
#  FC_FN_UNDERSCORE
#    "NO_UNDER" if nothing is appended, "UNDER" if
#    one underscore is appended, and "SECOND_UNDER"
#    if a function with an underscore has a second
#    appended.
#
# FC_FUNC_DEFAULT
#    The default mange mangling for Fortran functions
#    that do not contain an underscore.
#
# FC_FUNC__DEFAULT
#    The default mange mangling for Fortran functions
#    that do contain an underscore.
#
#  The Fortran 2003 name binding facilities and ISO_C_BINDING module
#  should be preferred over cpp macro trickery whenever possible.
#
function(fortran_mangling)

  if(NOT DEFINED FC_FN_CASE)

    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("FORTRAN_MANGLING: Testing name Mangling Schemes!\n")
    endif()

    find_file(_fcmakelists fmangle/ ${CMAKE_MODULE_PATH})
    if (NOT _fcmakelists)
      message(STATUS "Error, the directory fmangle could not be found so we can not determine Fortran name mangling!")
      return()
    endif()

    set(_fcmangledir ${PROJECT_BINARY_DIR}/CMakeFiles/CMakeTmp/fmangle)
    file(MAKE_DIRECTORY ${_fcmangledir})

    foreach(cdef LOWER UPPER)

      foreach(udef UNDER NO_UNDER SECOND_UNDER)

        if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          message("FORTRAN_MANGLING: Testing ${cdef} ${udef}\n\n")
        endif()

        set(_fcmangledir_case "${_fcmangledir}/${cdef}/${udef}")
        file(MAKE_DIRECTORY "${_fcmangledir}/${cdef}")
        file(MAKE_DIRECTORY ${_fcmangledir_case})

        set(COMMON_DEFS -DFC_FN_${cdef} -DFC_FN_${udef})
        set(C_FLAGS "${CMAKE_C_FLAGS} ${${PROJECT_NAME}_EXTRA_LINK_FLAGS}")
        set(F_FLAGS "${CMAKE_Fortran_FLAGS} ${${PROJECT_NAME}_EXTRA_LINK_FLAGS}")
        try_compile(_fcmngl ${_fcmangledir_case} ${_fcmakelists} fmangle
          CMAKE_FLAGS
            "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}"
            "-DCMAKE_C_FLAGS:STRING=${C_FLAGS}"
            "-DCMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}:STRING=${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}"
            "-DCMAKE_Fortran_FLAGS:STRING=${F_FLAGS}"
            "-DCMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}:STRING=${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}"
            "-DCOMMON_DEFS=${COMMON_DEFS}"
          OUTPUT_VARIABLE _fcmngl_output
          )
        if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
          message("${_fcmngl_output}\n\n")
        endif()

        if(_fcmngl)
          if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
            message("FORTRAN_MANGLING: Bingo!  ${cdef} ${udef} is the correct fortran name mangling!\n")
          endif()
          global_set(FC_FN_CASE ${cdef})
          global_set(FC_FN_UNDERSCORE ${udef})
          break()
        endif()

      endforeach()

      if(_fcmngl)
        break()
      endif()

    endforeach()

    if(_fcmngl)
      message(STATUS "Fortran name mangling: ${FC_FN_CASE} ${FC_FN_UNDERSCORE}")
    else()
      message(STATUS "Warning, cannot automatically determine Fortran mangling.")
    endif()

  endif()

  if (FC_FN_CASE STREQUAL LOWER)
    set(FC_NAME_NAME name)
  elseif (FC_FN_CASE STREQUAL UPPER)
    set(FC_NAME_NAME NAME)
  endif()

  if (FC_FN_UNDERSCORE)
    if(FC_FN_UNDERSCORE STREQUAL "UNDER")
      set(FC_FUNC_DEFAULT "(name,NAME) ${FC_NAME_NAME} ## _" CACHE INTERNAL "")
      set(FC_FUNC__DEFAULT "(name,NAME) ${FC_NAME_NAME} ## _" CACHE INTERNAL "")
    elseif(FC_FN_UNDERSCORE STREQUAL "SECOND_UNDER")
      set(FC_FUNC_DEFAULT "(name,NAME) ${FC_NAME_NAME} ## _" CACHE INTERNAL "")
      set(FC_FUNC__DEFAULT "(name,NAME) ${FC_NAME_NAME} ## __" CACHE INTERNAL "")
    else()
      set(FC_FUNC_DEFAULT "(name,NAME) ${FC_NAME_NAME}" CACHE INTERNAL "")
      set(FC_FUNC__DEFAULT "(name,NAME) ${FC_NAME_NAME}" CACHE INTERNAL "")
    endif()
  endif()

endfunction()


# 2008/10/26: rabartl: Below, these were macros that were also present in the
# file that I got off of the web but I did not need them so they are disabled
# (for now)

## - Guess if the Fortran compiler returns REAL in C doubles.
##  fortran_floatret()
##
## If the Fortran compiler follows the f2c convention of
## returning REALs in C doubles, routines like SLAMCH can
## cause very difficult to find stack corruptions.  This
## test is not perfect; it just checks if C->Fortran->C
## twice in a row does not crash.  If the tests do crash,
## FC_FN_FLOATRET is set to a true value, otherwise it
## is set to a false value.
##
## The REAL kinds in Fortran 2003's ISO_C_BINDING module should
## be used instead of this test whenever possible.
##
#macro(fortran_floatret)
#  if(NOT DEFINED FC_FN_FLOATRET)
#    # Find inputs
#    find_file(_fcindir floatret/ ${CMAKE_MODULE_PATH})
#    set(_fcdir ${PROJECT_BINARY_DIR}/CMakeFiles/CMakeTmp/floatret)
#    file(MAKE_DIRECTORY ${_fcdir})
#    try_compile(_fccrv ${_fcdir} ${_fcindir} floatret
#      CMAKE_FLAGS "-DFC_FN_DEFS:STRING=${FC_FN_DEFS}")
#    if(_fccrv)
#      execute_process(COMMAND ${_fcdir}/ctst
#        WORKING_DIRECTORY ${_fcdir}
#        RESULT_VARIABLE _fcrrv)
#    endif(_fccrv)
#    if(_fcrrv EQUAL 0)
#      set(_fc_fn_floatret 0)
#    else(_fcrrv EQUAL 0)
#      set(_fc_fn_floatret 1)
#    endif(_fcrrv EQUAL 0)
#    set(FC_FN_FLOATRET ${_fc_fn_floatret} CACHE BOOL
#      "Fortran returns REAL in double.")
#    message(STATUS "Fortran returns REAL in double: ${FC_FN_FLOATRET}")
#  endif(NOT DEFINED FC_FN_FLOATRET)
#endmacro()
#
#
## - Guess the convention for passing strings from C to Fortran.
##  fortran_stringarg()
##
## If string lengths are directly appended to each variable, e.g.
## CALL foo('bar', 1.0) becomes foo({'b','a','r'}, 3, 1.0), then
## FC_FN_STRINGARG is set to PAIRED.  If the lengths are appended
## to the call, e.g. foo({'b','a','r'}, 1.0, 3), FC_FN_STRINGARG
## is set to TRAILING.
##
## This macro does not currently check for older Cray and VMS
## conventions that require conversion functions.  It also assumes
## that the length is passed as the "natural" size type, C's size_t.
##
## The string kinds in Fortran 2003's ISO_C_BINDING module should
## be used instead of these conventions whenever possible.
##
#set(FC_FN_STRINGARG_TYPES "PAIRED" "TRAILING")
#macro(fortran_stringarg)
#  if(NOT DEFINED FC_FN_STRINGARG)
#    # Find inputs
#    find_file(_fcstrindir fstrings/ ${CMAKE_MODULE_PATH})
#    set(_fcstrdir ${PROJECT_BINARY_DIR}/CMakeFiles/CMakeTmp/fstrings)
#    foreach(argtype ${FC_FN_STRINGARG_TYPES})
#      file(MAKE_DIRECTORY ${_fcstrdir})
#      try_compile(_fcstrcrv ${_fcstrdir} ${_fcstrindir} fstrings
#        CMAKE_FLAGS "-DFC_FN_DEFS:STRING=${FC_FN_DEFS}" "-Dargtype:STRING=${argtype}")
#      if(_fcstrcrv)
#        execute_process(COMMAND ${_fcstrdir}/ccheck
#          WORKING_DIRECTORY ${_fcstrdir}
#          RESULT_VARIABLE _fcstrrrv)
#        if(_fcstrrrv EQUAL 0)
#          set(_fcstr ${argtype})
#          break()
#        endif(_fcstrrrv EQUAL 0)
#      endif(_fcstrcrv)
#      file(REMOVE_RECURSE ${_fcstrdir})
#    endforeach(argtype)
#    if(DEFINED _fcstr)
#      set(FC_FN_STRINGARG ${_fcstr} CACHE STRING
#        "How Fortran accepts string arguments.")
#      message(STATUS "Fortran string passing: ${FC_FN_STRINGARG}")
#    else(DEFINED _fcstr)
#      message(STATUS "Cannot determine Fortran string passing.")
#    endif(DEFINED _fcstr)
#  endif(NOT DEFINED FC_FN_STRINGARG)
#endmacro()

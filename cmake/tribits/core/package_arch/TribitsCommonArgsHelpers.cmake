# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

function(tribits_set_linker_language_from_arg  TARGET_NAME_IN  LINKER_LANGUAGE_IN)

  if( TRIBITS_SET_LINKER_LANGUAGE_FROM_ARG_DEBUG_DUMP)
    message("tribits_set_linker_language_from_arg( '${TARGET_NAME_IN}' '${LINKER_LANGUAGE_IN}' )")
  endif()

  if(LINKER_LANGUAGE_IN)
    set(LINKER_LANGUAGE ${LINKER_LANGUAGE_IN})
  elseif (${PROJECT_NAME}_ENABLE_CXX)
    set(LINKER_LANGUAGE CXX)
  elseif(${PROJECT_NAME}_ENABLE_C)
    set(LINKER_LANGUAGE C)
  else()
    set(LINKER_LANGUAGE)
  endif()

  if (LINKER_LANGUAGE)

    if (${PROJECT_NAME}_VERBOSE_CONFIGURE
        OR TRIBITS_SET_LINKER_LANGUAGE_FROM_ARG_DEBUG_DUMP
      )
      message("-- Setting linker language for target '${TARGET_NAME_IN}' to '${LINKER_LANGUAGE}'")
    endif()

    set_property(
      TARGET ${TARGET_NAME_IN}
      APPEND PROPERTY LINKER_LANGUAGE ${LINKER_LANGUAGE}
      )

  endif()

endfunction()



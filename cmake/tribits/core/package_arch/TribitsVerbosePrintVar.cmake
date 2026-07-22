# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


if (${PROJECT_NAME}_VERBOSE_CONFIGURE)


# @FUNCTION: tribits_verbose_print_var()
#
# print a variable giving its name then value if
# ``${PROJECT_NAME}_VERBOSE_CONFIGURE=TRUE``.
#
# Usage::
#
#   tribits_verbose_print_var(<varName>)
#
# This prints::
#
#   message("-- " "${VARIBLE_NAME}='${${VARIBLE_NAME}}'")
#
# The variable ``<varName>`` can be defined or undefined or empty.  This uses
# an explicit "-- " line prefix so that it prints nice even on Windows CMake.
#
function(tribits_verbose_print_var VARIBLE_NAME)
  message("-- " "${VARIBLE_NAME}='${${VARIBLE_NAME}}'")
endfunction()


else() # ${PROJECT_NAME}_VERBOSE_CONFIGURE


function(tribits_verbose_print_var VARIBLE_NAME)
endfunction()


endif() # ${PROJECT_NAME}_VERBOSE_CONFIGURE

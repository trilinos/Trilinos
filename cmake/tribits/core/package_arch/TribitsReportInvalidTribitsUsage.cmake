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

IF (__TribitsReportInvalidTribitsUsage_INCLUDED__)
  RETURN()
ELSE()
  SET(__TribitsReportInvalidTribitsUsage_INCLUDED__ TRUE)
ENDIF()

INCLUDE(MessageWrapper)

# Called to report incorrect usage
#
# Usage:
#
#   TRIBITS_REPORT_INVALID_TRIBITS_USAGE("<err_msg0>" "<err_msg1>" ...)
#
# Depending on the value of ${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE, this
# function will:
#
#   * FATAL_ERROR: Calls MESSAGE(FATAL_ERROR "<error_message>")
#   * SEND_ERROR: Calls MESSAGE(SEND_ERROR "<error_message>")
#   * WARNING: Calls MESSAGE(WARNING "<error_message>")
#   * IGNORE: Does not call MESSAGE() at all and is silent
#
# If '${${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE}' is empty on call, then
# `FATAL_ERROR` will be used.
#
FUNCTION(TRIBITS_REPORT_INVALID_TRIBITS_USAGE)
  IF ("${${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE}" STREQUAL "")
    SET(${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE FATAL_ERROR)
  ENDIF()
  SET(PRINT_ERR_MSG)
  IF (${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE STREQUAL "FATAL_ERROR")
    SET(PRINT_ERR_MSG FATAL_ERROR)
  ELSEIF (${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE STREQUAL "SEND_ERROR")
    SET(PRINT_ERR_MSG SEND_ERROR)
  ELSEIF (${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE STREQUAL "WARNING")
    SET(PRINT_ERR_MSG WARNING)
  ELSEIF (${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE STREQUAL "IGNORE")
    SET(PRINT_ERR_MSG)
  ELSE()
    MESSAGE_WRAPPER(FATAL_ERROR "Error, invalid value for"
      " ${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE ="
      " '${${PROJECT_NAME}_ASSERT_CORRECT_TRIBITS_USAGE}'!"
      "  Value values include 'FATAL_ERROR', 'SEND_ERROR', 'WARNING', and 'IGNORE'!")
  ENDIF()
  IF (PRINT_ERR_MSG)
    MESSAGE_WRAPPER(${PRINT_ERR_MSG} ${ARGN})
  ENDIF()
ENDFUNCTION()

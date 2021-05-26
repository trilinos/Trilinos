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


# @FUNCTION: TRIBITS_PRINT_INITIAL_DEPENDENCY_INFO()
#
# Usage::
#
#   TRIBITS_PRINT_INITIAL_DEPENDENCY_INFO()
#
# Function that prints whatever initial dependency information that is
# available that is requested by the user after the initial construction of
# the package dependency graph but **before** the call of
# `TRIBITS_ADJUST_PACKAGE_ENABLES()`_.
#
# Calls:
#
# * `TRIBITS_PRINT_TENTATIVELY_ENABLED_TPLS()`_
# * `TRIBITS_DUMP_PACKAGE_DEPENDENCIES_INFO()`_
#
FUNCTION(TRIBITS_PRINT_INITIAL_DEPENDENCY_INFO)
  TRIBITS_PRINT_TENTATIVELY_ENABLED_TPLS()
  TRIBITS_DUMP_PACKAGE_DEPENDENCIES_INFO()
ENDFUNCTION()


# @FUNCTION: TRIBITS_PRINT_TENTATIVELY_ENABLED_TPLS()
#
# Usage::
#
#   TRIBITS_PRINT_TENTATIVELY_ENABLED_TPLS()
#
# Function that print the set of tentatively enabled TPLs.
#
# Does **not** modify any state!
#
FUNCTION(TRIBITS_PRINT_TENTATIVELY_ENABLED_TPLS)
  FOREACH(TPL ${${PROJECT_NAME}_TPLS})
    IF (TPL_TENTATIVE_ENABLE_${TPL})
      MESSAGE("-- Tentatively enabling TPL '${TPL}'")
      #PRINT_VAR(TPL_ENABLE_${TPL})
    ENDIF()
  ENDFOREACH()
ENDFUNCTION()


# @FUNCTION: TRIBITS_DUMP_PACKAGE_DEPENDENCIES_INFO()
#
# Usage:
#
#  TRIBITS_DUMP_PACKAGE_DEPENDENCIES_INFO()
#
# Function that dumps (prints to STDOUT) the package dependency info if
# ``${PROJECT_NAME}_DUMP_PACKAGE_DEPENDENCIES==TRUE``.
#
# Does **not** modify state!
#
FUNCTION(TRIBITS_DUMP_PACKAGE_DEPENDENCIES_INFO)

  ADVANCED_OPTION(${PROJECT_NAME}_DUMP_PACKAGE_DEPENDENCIES
    "Dump the package dependency information."
    "${${PROJECT_NAME}_VERBOSE_CONFIGURE}" )

  ADVANCED_OPTION(${PROJECT_NAME}_DUMP_FORWARD_PACKAGE_DEPENDENCIES
    "Dump the package forward dependency information."
    "${${PROJECT_NAME}_VERBOSE_CONFIGURE}" )

  IF (${PROJECT_NAME}_DUMP_PACKAGE_DEPENDENCIES)
    MESSAGE("")
    MESSAGE("Printing package dependencies ...")
    MESSAGE("")
    PRINT_NONEMPTY_VAR_WITH_SPACES(${PROJECT_NAME}_PACKAGES  DUMMY_OUT)
    MESSAGE("")
    PRINT_NONEMPTY_VAR_WITH_SPACES(${PROJECT_NAME}_SE_PACKAGES  DUMMY_OUT)
    MESSAGE("")
    FOREACH(TRIBITS_PACKAGE ${${PROJECT_NAME}_SE_PACKAGES})
      TRIBITS_PRINT_PACKAGE_DEPENDENCIES(${TRIBITS_PACKAGE})
      MESSAGE("")
    ENDFOREACH()
  ENDIF()

ENDFUNCTION()

#  LocalWords:  TRIBITS

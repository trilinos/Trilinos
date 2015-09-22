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


# Standard TriBITS Includes
INCLUDE(TribitsTplFindIncludeDirsAndLibraries)
INCLUDE(TribitsGeneralMacros)

# Standard TriBITS utilities includes
INCLUDE(AppendStringVar)


#
# @FUNCTION: TRIBITS_PROCESS_ENABLED_TPL()
#
# Processs an enabled TPL's FindTPL${TPL_NAME}.cmake module.
#
FUNCTION(TRIBITS_PROCESS_ENABLED_TPL  TPL_NAME)

  # Setup the processing string
  SET(PROCESSING_MSG_STRING "Processing enabled TPL: ${TPL_NAME} (")
  IF (TPL_${TPL_NAME}_ENABLING_PKG)
    APPEND_STRING_VAR(PROCESSING_MSG_STRING
      "enabled by ${TPL_${TPL_NAME}_ENABLING_PKG}," )
  ELSE()
    APPEND_STRING_VAR(PROCESSING_MSG_STRING
      "enabled explicitly," )
  ENDIF()
    APPEND_STRING_VAR(PROCESSING_MSG_STRING 
      " disable with -DTPL_ENABLE_${TPL_NAME}=OFF)" )

  # Print the processing header
  MESSAGE("${PROCESSING_MSG_STRING}")

  IF (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)

    # Locate the FindTPL${TPL_NAME}.cmake module
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      PRINT_VAR(${TPL_NAME}_FINDMOD)
      PRINT_VAR(${TPL_NAME}_TPLS_LIST_FILE)
    ENDIF()
    IF (IS_ABSOLUTE ${${TPL_NAME}_FINDMOD})
      #MESSAGE("${${TPL_NAME}_FINDMOD} is absolute!")
      SET(CURRENT_TPL_PATH "${${TPL_NAME}_FINDMOD}")
    ELSE()
      #MESSAGE("${${TPL_NAME}_FINDMOD} is *NOT* absolute!")
      SET(CURRENT_TPL_PATH "${PROJECT_SOURCE_DIR}/${${TPL_NAME}_FINDMOD}")
    ENDIF()
    #PRINT_VAR(CURRENT_TPL_PATH)

    # Process the FindTPL${TPL_NAME}.cmake module
    TRIBITS_TRACE_FILE_PROCESSING(TPL  INCLUDE  "${CURRENT_TPL_PATH}")
    INCLUDE("${CURRENT_TPL_PATH}")

    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      PRINT_VAR(TPL_${TPL_NAME}_NOT_FOUND)
    ENDIF()

    # Address failed find of the TPL
    IF (TPL_${TPL_NAME}_NOT_FOUND AND NOT TPL_TENTATIVE_ENABLE_${TPL_NAME})
      MESSAGE(
        "-- NOTE: The find module file for this failed TPL '${TPL_NAME}' is:\n"
        "     ${CURRENT_TPL_PATH}\n"
        "   which is pointed to in the file:\n"
        "     ${${TPL_NAME}_TPLS_LIST_FILE}\n"
        )
      IF (TPL_${TPL_NAME}_ENABLING_PKG)
        MESSAGE(
          "TIP: One way to get past the configure failure for the\n"
          "TPL '${TPL_NAME}' is to simply disable it with:\n"
          "  -DTPL_ENABLE_${TPL_NAME}=OFF\n"
          "which will disable it and will recursively disable all of the\n"
          "downstream packages that have required dependencies on it, including\n"
          "the package '${TPL_${TPL_NAME}_ENABLING_PKG}' which triggered its enable.\n"
          "When you reconfigure, just grep the cmake stdout for '${TPL_NAME}'\n"
          "and then follow the disables that occur as a result to see what impact\n"
          "this TPL disable has on the configuration of ${PROJECT_NAME}.\n"
          )
      ELSE()
        MESSAGE(
          "TIP: Even though the TPL '${TPL_NAME}' was explicitly enabled in input,\n"
          "it can be disabled with:\n"
          "  -DTPL_ENABLE_${TPL_NAME}=OFF\n"
          "which will disable it and will recursively disable all of the\n"
          "downstream packages that have required dependencies on it.\n"
          "When you reconfigure, just grep the cmake stdout for '${TPL_NAME}'\n"
          "and then follow the disables that occur as a result to see what impact\n"
          "this TPL disable has on the configuration of ${PROJECT_NAME}.\n"
          )
      ENDIF()
      MESSAGE(FATAL_ERROR
        "ERROR: TPL_${TPL_NAME}_NOT_FOUND=${TPL_${TPL_NAME}_NOT_FOUND}, aborting!")
    ENDIF()

    # Assert that the TPL correctly defined all of these varaibles.
    ASSERT_DEFINED(TPL_${TPL_NAME}_INCLUDE_DIRS)
    ASSERT_DEFINED(TPL_${TPL_NAME}_LIBRARIES)
    ASSERT_DEFINED(TPL_${TPL_NAME}_LIBRARY_DIRS)
    # ToDo: Make TPL_${TPL_NAME}_LIBRARY_DIRS go away.  It is not needed for
    # anything.

  ENDIF()

ENDFUNCTION()

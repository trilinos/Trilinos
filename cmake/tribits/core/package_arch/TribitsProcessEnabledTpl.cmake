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
include(TribitsExternalPackageFindTplHelpers)
include(TribitsExternalPackageWriteConfigFile)
include(TribitsTplFindIncludeDirsAndLibraries)
include(TribitsGeneralMacros)

# Standard TriBITS utilities includes
include(AppendStringVar)
include(TribitsStandardizePaths)


# @FUNCTION: tribits_process_enabled_tpl()
#
# Process an enabled TPL's FindTPL${TPL_NAME}.cmake module.
#
function(tribits_process_enabled_tpl  TPL_NAME)

  # Setup the processing string
  set(PROCESSING_MSG_STRING "Processing enabled external package/TPL: ${TPL_NAME} (")
  if (${TPL_NAME}_ENABLING_PKG)
    string(APPEND PROCESSING_MSG_STRING
      "enabled by ${${TPL_NAME}_ENABLING_PKG}," )
  else()
    string(APPEND PROCESSING_MSG_STRING
      "enabled explicitly," )
  endif()
    string(APPEND PROCESSING_MSG_STRING 
      " disable with -DTPL_ENABLE_${TPL_NAME}=OFF)" )

  # Print the processing header
  message("${PROCESSING_MSG_STRING}")

  if (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)

    # Locate the FindTPL${TPL_NAME}.cmake module
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      print_var(${TPL_NAME}_FINDMOD)
    endif()
    if (${TPL_NAME}_FINDMOD STREQUAL "TRIBITS_PKG")
      set(TPL_${TPL_NAME}_PARTS_ALREADY_SET FALSE)  # ToDo: Take out?
      if (NOT TPL_${TPL_NAME}_PARTS_ALREADY_SET)
        find_package(${TPL_NAME} CONFIG REQUIRED)
        global_set(TPL_${TPL_NAME}_LIBRARIES
          "${${TPL_NAME}_LIBRARIES}" "${${TPL_NAME}_TPL_LIBRARIES}")
        global_set(TPL_${TPL_NAME}_PARTS_ALREADY_SET TRUE)
      endif()
      return()
    elseif (IS_ABSOLUTE ${${TPL_NAME}_FINDMOD})
      #message("${${TPL_NAME}_FINDMOD} is absolute!")
      set(CURRENT_TPL_PATH "${${TPL_NAME}_FINDMOD}")
    else()
      #message("${${TPL_NAME}_FINDMOD} is *NOT* absolute!")
      set(CURRENT_TPL_PATH "${PROJECT_SOURCE_DIR}/${${TPL_NAME}_FINDMOD}")
    endif()
    #print_var(CURRENT_TPL_PATH)

    # Process the FindTPL${TPL_NAME}.cmake module
    tribits_trace_file_processing(TPL  INCLUDE  "${CURRENT_TPL_PATH}")
    set(TRIBITS_FINDING_RAW_${TPL_NAME}_PACKAGE_FIRST TRUE)
    include("${CURRENT_TPL_PATH}")
    unset(TRIBITS_FINDING_RAW_${TPL_NAME}_PACKAGE_FIRST)
    # NOTE: Above, setting TRIBITS_FINDING_RAW_${TPL_NAME}_PACKAGE_FIRST=TRUE
    # triggers special logic in the TriBITS-created
    # ${TPL_NAME}ConfigVersion.cmake file to set
    # PACKAGE_VERSION_COMPATIBLE=FALSE and result in find_package(${TPL_NAME})
    # that may be called inside of ${TPL_NAME}_FINDMOD to not find a
    # TriBITS-generated ${TPL_NAME}Config.cmake file.  This allows
    # find_package(${TPL_NAME}) to usae a proper non-TriBITS
    # Find${TPL_NAME}.cmake module or find a non-TriBITS
    # ${TPL_NAME}Config.cmake module.

    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      print_var(TPL_${TPL_NAME}_NOT_FOUND)
    endif()

    # Address failed find of the TPL
    if (TPL_${TPL_NAME}_NOT_FOUND AND NOT TPL_TENTATIVE_ENABLE_${TPL_NAME})
      message(
        "-- NOTE: The find module file for this failed TPL '${TPL_NAME}' is:\n"
        "     ${CURRENT_TPL_PATH}\n"
        "   which is pointed to in the file:\n"
        "     ${${TPL_NAME}_TPLS_LIST_FILE}\n"
        )
      if (${TPL_NAME}_ENABLING_PKG)
        message(
          "TIP: One way to get past the configure failure for the\n"
          "TPL '${TPL_NAME}' is to simply disable it with:\n"
          "  -DTPL_ENABLE_${TPL_NAME}=OFF\n"
          "which will disable it and will recursively disable all of the\n"
          "downstream packages that have required dependencies on it, including\n"
          "the package '${${TPL_NAME}_ENABLING_PKG}' which triggered its enable.\n"
          "When you reconfigure, just grep the cmake stdout for '${TPL_NAME}'\n"
          "and then follow the disables that occur as a result to see what impact\n"
          "this TPL disable has on the configuration of ${PROJECT_NAME}.\n"
          )
      else()
        message(
          "TIP: Even though the TPL '${TPL_NAME}' was explicitly enabled in input,\n"
          "it can be disabled with:\n"
          "  -DTPL_ENABLE_${TPL_NAME}=OFF\n"
          "which will disable it and will recursively disable all of the\n"
          "downstream packages that have required dependencies on it.\n"
          "When you reconfigure, just grep the cmake stdout for '${TPL_NAME}'\n"
          "and then follow the disables that occur as a result to see what impact\n"
          "this TPL disable has on the configuration of ${PROJECT_NAME}.\n"
          )
      endif()
      message(FATAL_ERROR
        "ERROR: TPL_${TPL_NAME}_NOT_FOUND=${TPL_${TPL_NAME}_NOT_FOUND}, aborting!")
    endif()

    # Generate the <tplName>ConfigVersion.cmake file if it has not been
    # created yet and add install targets for <tplName>Config[Version].cmake
    set(buildDirExternalPkgsDir
      "${${PROJECT_NAME}_BINARY_DIR}/${${PROJECT_NAME}_BUILD_DIR_EXTERNAL_PKGS_DIR}")
    set(tplConfigFile
      "${buildDirExternalPkgsDir}/${TPL_NAME}/${TPL_NAME}Config.cmake")
    set(tplConfigVersionFile
      "${buildDirExternalPkgsDir}/${TPL_NAME}/${TPL_NAME}ConfigVersion.cmake")
    tribits_extpkg_write_config_version_file(${TPL_NAME}
      "${tplConfigVersionFile}")
    tribits_extpkg_install_config_file(${TPL_NAME} "${tplConfigFile}")
    tribits_extpkg_install_config_version_file(${TPL_NAME}
      "${tplConfigVersionFile}")

  endif()

endfunction()

# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# Standard TriBITS Includes
include(TribitsExternalPackageWithImportedTargetsFindTplModuleHelpers)
include(TribitsExternalPackageWriteConfigFile)
include(TribitsTplFindIncludeDirsAndLibraries)
include(TribitsGeneralMacros)
include(TribitsConfigureTiming)

# Standard TriBITS utilities includes
include(AppendStringVar)
include(TribitsStandardizePaths)
include(TribitsCreateReverseList)


# @MACRO: tribits_process_enabled_tpls()
#
# Gather information and targets from enabled TPLs
#
# For more info, see `Processing of external packages/TPLs and
# TriBITS-compliant external packages`_.
#
macro(tribits_process_enabled_tpls)

  tribits_config_code_start_timer(CONFIGURE_TPLS_TIME_START_SECONDS)

  tribits_filter_package_list_from_var(${PROJECT_NAME}_DEFINED_TOPLEVEL_PACKAGES
    EXTERNAL  ON  NONEMPTY  ${PROJECT_NAME}_enabledExternalTopLevelPackages)

  tribits_project_has_tribits_compliant_external_packages(
    ${PROJECT_NAME}_enabledExternalTopLevelPackages
    projectHasTribitsCompliantExternalPackages )

  if (projectHasTribitsCompliantExternalPackages)
    message("")
    message("Getting information for all enabled TriBITS-compliant"
      " or upstream external packages/TPLs in reverse order ...")
    message("")

    tribits_create_reverse_list(${PROJECT_NAME}_enabledExternalTopLevelPackages
      ${PROJECT_NAME}_reverseEnabledExternalTopLevelPackages)

    foreach(TPL_NAME  IN LISTS  ${PROJECT_NAME}_reverseEnabledExternalTopLevelPackages)
      if (${TPL_NAME}_IS_TRIBITS_COMPLIANT
          OR ${TPL_NAME}_PROCESSED_BY_DOWNSTREAM_TRIBITS_EXTERNAL_PACKAGE
        )
        tribits_process_enabled_tribits_compliant_or_upstream_tpl(${TPL_NAME})
      endif()
    endforeach()

    set(remainingTplsTextStr " remaining")

  else()

    set(remainingTplsTextStr "")

  endif()

  message("")
  message("Getting information for all${remainingTplsTextStr} enabled external packages/TPLs ...")
  message("")

  foreach(TPL_NAME  IN LISTS  ${PROJECT_NAME}_enabledExternalTopLevelPackages)
    if ((NOT ${TPL_NAME}_IS_TRIBITS_COMPLIANT)
        AND (NOT ${TPL_NAME}_PROCESSED_BY_DOWNSTREAM_TRIBITS_EXTERNAL_PACKAGE)
      )
      tribits_process_enabled_standard_tpl(${TPL_NAME})
    endif()
  endforeach()

  tribits_config_code_stop_timer(CONFIGURE_TPLS_TIME_START_SECONDS
    "\nTotal time to configure enabled external packages/TPLs")

endmacro()


macro(tribits_process_enabled_tribits_compliant_or_upstream_tpl  TPL_NAME)

  tribits_get_enabled_tpl_processing_string(${TPL_NAME}  tplProcessingString)
  message("${tplProcessingString}")

  if (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)
    if ( (NOT TARGET ${TPL_NAME}::all_libs)  AND  ${TPL_NAME}_IS_TRIBITS_COMPLIANT )
      tribits_process_enabled_tribits_compliant_tpl(${TPL_NAME})
      set(${TPL_NAME}_PROCESSED_BY_DOWNSTREAM_TRIBITS_EXTERNAL_PACKAGE  FALSE)
    elseif (TARGET ${TPL_NAME}::all_libs)
      message("-- "
        "The external package/TPL ${TPL_NAME} was defined by a downstream"
        " TriBITS-compliant external package already processed")
    elseif (${TPL_NAME}_FINDMOD AND (NOT ${TPL_NAME}_FINDMOD STREQUAL "TRIBITS_PKG"))
      message("-- "
        "The external package/TPL ${TPL_NAME} was *NOT* defined by a downstream"
        " TriBITS-compliant external package and must be found again in below loop")
      set(${TPL_NAME}_PROCESSED_BY_DOWNSTREAM_TRIBITS_EXTERNAL_PACKAGE  FALSE)
    else()
      message(FATAL_ERROR
        "Error, the external package/TPL ${TPL_NAME} was *NOT* defined by a downstream"
        " TriBITS-compliant external package and has not find module!")
    endif()
  endif()

endmacro()
# NOTE: Above, handles the case where an upstream external package/TPL should
# have been defined a downstream external package that was already processed
# but it was not defined (because the downstream packages was not a fully
# TriBITS-compliant external package).  For a TriBITS-compliant external
# package/TPL that should have been defined by a downstream TriBITS-compliant
# an external package/TPL, the first if-statement above takes care of that
# case by calling find_package(${TPL_NAME}) (because ${TPL_NAME}::all_libs
# will not be defined).  However, if the upstream external package/TPL is
# *NOT* TriBITS-compliant, then it may be a legacy TriBITS TPL which means
# that it must be processed in ascending order in order to build the
# downstream TriBITS TPLs correctly.


# @MACRO: tribits_process_enabled_standard_tpl()
#
# Process an enabled TPL's FindTPL${TPL_NAME}.cmake module.
#
macro(tribits_process_enabled_standard_tpl  TPL_NAME)

  tribits_get_enabled_tpl_processing_string(${TPL_NAME}  tplProcessingString)
  message("${tplProcessingString}")

  if (NOT ${PROJECT_NAME}_TRACE_DEPENDENCY_HANDLING_ONLY)

    # Locate the FindTPL${TPL_NAME}.cmake module
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      print_var(${TPL_NAME}_FINDMOD)
    endif()

    tribits_process_enabled_tribits_find_tpl_mod_file(${TPL_NAME})
    tribits_address_failed_tpl_find(${TPL_NAME})
    tribits_generate_tpl_version_file_and_add_package_config_install_targets(
      ${TPL_NAME})

  endif()

endmacro()


# Get external package/TPL processing string
#
function(tribits_get_enabled_tpl_processing_string  TPL_NAME  tplProcessingStringOut)
  set(tplProcessingString "Processing enabled external package/TPL: ${TPL_NAME} (")
  if (${TPL_NAME}_ENABLING_PKG)
    string(APPEND tplProcessingString "enabled by ${${TPL_NAME}_ENABLING_PKG}," )
  else()
    string(APPEND tplProcessingString "enabled explicitly," )
  endif()
  string(APPEND tplProcessingString " disable with -DTPL_ENABLE_${TPL_NAME}=OFF)" )
  set(${tplProcessingStringOut} "${tplProcessingString}" PARENT_SCOPE)
endfunction()


# Process an enabled TPL defined using a TriBITS-compliant external package
# <tplName>Config.cmake file
#
macro(tribits_process_enabled_tribits_compliant_tpl  TPL_NAME)
  message("-- "
    "Calling find_package(${TPL_NAME}) for TriBITS-compliant external package")
  find_package(${TPL_NAME} CONFIG REQUIRED)
  if (${TPL_NAME}_DIR)
    message("-- " "Found ${TPL_NAME}_DIR='${${TPL_NAME}_DIR}'")
  else()
    message(FATAL_ERROR
      "ERROR! Failed to find TriBITS-compliant external package ${TPL_NAME}!")
  endif()
endmacro()


# Process an enabled TPL defined using a FindTPL<tplName>.cmake module
#
macro(tribits_process_enabled_tribits_find_tpl_mod_file  TPL_NAME)

    if (IS_ABSOLUTE ${${TPL_NAME}_FINDMOD})
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

endmacro()


function(tribits_address_failed_tpl_find  TPL_NAME)
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
endfunction()


# Generate the <tplName>ConfigVersion.cmake file for a TriBITS TPL and install
# the already generated <tplName>Config.cmake file
#
function(tribits_generate_tpl_version_file_and_add_package_config_install_targets
    TPL_NAME
  )
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
endfunction()


function(tribits_project_has_tribits_compliant_external_packages
    enabledExternalTopLevelPackagesListName
    projectHasTribitsCompliantExternalPackagesOut
  )

  set(projectHasTribitsCompliantExternalPackages FALSE)

  foreach(TPL_NAME  IN LISTS  ${enabledExternalTopLevelPackagesListName})
    if (${TPL_NAME}_IS_TRIBITS_COMPLIANT
        OR ${TPL_NAME}_PROCESSED_BY_DOWNSTREAM_TRIBITS_EXTERNAL_PACKAGE
      )
      set(projectHasTribitsCompliantExternalPackages TRUE)
      break()
    endif()
  endforeach()

  set(${projectHasTribitsCompliantExternalPackagesOut}
    ${projectHasTribitsCompliantExternalPackages} PARENT_SCOPE)

endfunction()

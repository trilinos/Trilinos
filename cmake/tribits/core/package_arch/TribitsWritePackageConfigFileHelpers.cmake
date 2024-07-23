# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @FUNCTION: tribits_extpkg_append_tribits_compliant_package_config_vars_str()
#
# Append the standard TriBITS-compliant external package variables
#
function(tribits_extpkg_append_tribits_compliant_package_config_vars_str  packageName
    packageConfigCodeStrVarInOut
  )
  set(packageConfigCodeStr "${${packageConfigCodeStrVarInOut}}")
  string(APPEND packageConfigCodeStr
    "\n# Standard TriBITS-compliant external package variables\n"
    "set(${packageName}_IS_TRIBITS_COMPLIANT TRUE)\n"
    "set(${packageName}_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE \"\${CMAKE_CURRENT_LIST_FILE}\")\n"
    "set(${packageName}_TRIBITS_COMPLIANT_PACKAGE_CONFIG_FILE_DIR \"\${CMAKE_CURRENT_LIST_DIR}\")\n"
    )
  set(${packageConfigCodeStrVarInOut} "${packageConfigCodeStr}" PARENT_SCOPE)
endfunction()


function(tribits_extpkg_write_package_config_file_from_str  tplName  configFileStr)
  set(buildDirExternalPkgsDir
    "${${PROJECT_NAME}_BINARY_DIR}/${${PROJECT_NAME}_BUILD_DIR_EXTERNAL_PKGS_DIR}")
  set(tplConfigFile
    "${buildDirExternalPkgsDir}/${tplName}/${tplName}Config.cmake")
  file(WRITE "${tplConfigFile}" "${configFileStr}")
endfunction()
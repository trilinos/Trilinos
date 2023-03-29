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
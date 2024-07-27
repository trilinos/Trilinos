# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

# Define the TriBITS minimum required CMake version

set(TRIBITS_CMAKE_MINIMUM_REQUIRED 3.23.0)

macro(tribits_asesrt_minimum_cmake_version)

  if (CMAKE_VERSION VERSION_LESS ${TRIBITS_CMAKE_MINIMUM_REQUIRED})
    message(FATAL_ERROR "Error, TriBiTS must have version"
        " ${TRIBITS_CMAKE_MINIMUM_REQUIRED} or higher!")
  endif()
  
endmacro()

# File names for TriBITS system

set(${PROJECT_NAME}_PACKAGES_FILE_NAME PackagesList.cmake)

set(${PROJECT_NAME}_TPLS_FILE_NAME TPLsList.cmake)

set(${PROJECT_NAME}_EXTRA_EXTERNAL_REPOS_FILE_NAME ExtraRepositoriesList.cmake)

set(${PROJECT_NAME}_EXTRA_PACKAGES_FILE_NAME PackagesList.cmake)

set(${PROJECT_NAME}_REPO_VERSION_FILE_NAME ${PROJECT_NAME}RepoVersion.txt)

set(${PROJECT_NAME}_EXTRA_TPLS_FILE_NAME TPLsList.cmake)

# Directories relative to the TriBITS base directory

set(TRIBITS_PYTHON_SCRIPTS_DIR "python_utils")

set(TRIBITS_CI_SUPPORT_DIR "ci_support")

set(TRIBITS_CTEST_DRIVER_DIR "ctest_driver")

set(TRIBITS_CMAKE_UTILS_DIR "core/utils")

set(TRIBITS_CMAKE_PACKAGE_ARCH_DIR "core/package_arch")

set(TRIBITS_CMAKE_INSTALLATION_FILES_DIR "core/installation")

# Files and directories related to the specific project

set(${PROJECT_NAME}_PACKAGE_DEPS_XML_FILE_NAME ${PROJECT_NAME}PackageDependencies.xml)

set(${PROJECT_NAME}_CDASH_SUBPROJECT_DEPS_XML_FILE_NAME CDashSubprojectDependencies.xml)

set(${PROJECT_NAME}_PACKAGE_DEPS_TABLE_HTML_FILE_NAME ${PROJECT_NAME}PackageDependenciesTable.html)

set(${PROJECT_NAME}_PACKAGE_DEPS_FILES_DIR "cmake/dependencies")

set(${PROJECT_NAME}_BUILD_DIR_EXTERNAL_PKGS_DIR "external_packages")

set(${PROJECT_NAME}_BUILD_DIR_CMAKE_PKGS_DIR "cmake_packages")

# Other stuff

if(${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Windows")
  #Apparently FIND_PROGRAM looks for an exact match of the file name.
  #So even though "git clone ..." is valid to use on windows we need to give the
  #full name of the command we want to run.
  set(GIT_NAME git.cmd)
else()
  set(GIT_NAME git)
endif()

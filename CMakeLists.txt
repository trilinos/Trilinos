#
# A) Define your project name and set up major project options
#

# To be safe, define your minimum CMake version.  This may be newer than the
# min required by TriBITS.
CMAKE_MINIMUM_REQUIRED(VERSION 3.23.0 FATAL_ERROR)

# Must set the project name as a variable at very beginning before including anything else
# We set the project name in a separate file so CTest scripts can use it.
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/ProjectName.cmake)

# CMake requires that you declare the CMake project in the top-level file and
# not in an include file :-(
PROJECT(${PROJECT_NAME} NONE)

# Set up to use ccache
include("${CMAKE_CURRENT_LIST_DIR}/cmake/UseCCache.cmake")

set(DEPRECATED_PACKAGES Amesos AztecOO Epetra EpetraExt Ifpack Intrepid Isorropia ML NewPackage Pliris PyTrilinos ShyLU_DDCore ThyraEpetraAdapters ThyraEpetraExtAdapters Triutils)

# Set an env so we know we are in configure
set(ENV{CMAKE_IS_IN_CONFIGURE_MODE} 1)

# Don't define TriBITS override of include_directories()
set(TRIBITS_HIDE_DEPRECATED_INCLUDE_DIRECTORIES_OVERRIDE TRUE)

#
# B) Pull in the TriBITS system and execute
#

INCLUDE(${CMAKE_CURRENT_LIST_DIR}/cmake/tribits/TriBITS.cmake)

# Make Trilinos create <Package>Config.cmake files by default
SET(${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES_DEFAULT ON)
# Make Trilinos set up CPack support by default
SET(${PROJECT_NAME}_ENABLE_CPACK_PACKAGING_DEFAULT ON)
# Don't allow disabled subpackages to be excluded from tarball
SET(${PROJECT_NAME}_EXCLUDE_DISABLED_SUBPACKAGES_FROM_DISTRIBUTION_DEFAULT FALSE)


SET(Trilinos_USE_GNUINSTALLDIRS_DEFAULT ON)

SET(Trilinos_MUST_FIND_ALL_TPL_LIBS_DEFAULT TRUE)

# Some CMake and TriBiTS tweaks just for Trilinos
include(TrilinosTweaks)

# Do all of the processing for this Tribits project
TRIBITS_PROJECT()

INSTALL_BUILD_STATS_SCRIPTS()

# Install TriBITS so that other projects can use it
include(SetupTribitsInstall)

IF(${PROJECT_NAME}_ENABLE_YouCompleteMe)
  INCLUDE(CodeCompletion)
ENDIF()

set(enabled_deprecated_packages "")
FOREACH(package ${DEPRECATED_PACKAGES})
  IF(Trilinos_ENABLE_${package})
    set(enabled_deprecated_packages ${package} ${enabled_deprecated_packages})
  ENDIF()
ENDFOREACH()

IF(enabled_deprecated_packages)
  message("********************************************************************")
  message("********  WARNING WARNING WARNING WARNING WARNING WARNING  *********")
  message("********************************************************************")
  message("* Several Trilinos packages are scheduled for deprecation in 2025.")
  message("* The following deprecated packages are enabled (either explicitly")
  message("* or implicitly due to a dependency)")
  message("*")
  string(REPLACE ";" "\n*    " p "${enabled_deprecated_packages}")
  message("*    ${p}")
  message("*")
  message("********************************************************************")
  message("")
  message(WARNING "Deprecated packages that are enabled: ${enabled_deprecated_packages}")
  IF(Trilinos_HIDE_DEPRECATED_CODE)
    message(FATAL_ERROR "Trilinos cannot be built with deprecated code hidden if any deprecated packages are enabled!")
  ENDIF()
ENDIF()

message(STATUS "If publishing results using Trilinos, please cite us: https://trilinos.github.io/cite.html")

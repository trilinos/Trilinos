# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


################################################################################
#
# Module TribitsExternalPackageWithImportedTargetsFindTplModuleHelpers.cmake
#
# Contains functions for implementing FindTPL<tplName>.cmake files for
# external packages using find_package(<externalPkg>) that producing modern
# IMPORTED targets.
#
# NOTE: The acronym 'extpkgwit' stands for "External Package With Imported
# Targets".
# 
################################################################################


include(TribitsExternalPackageWriteConfigFile)



# @FUNCTION: tribits_extpkg_create_imported_all_libs_target_and_config_file()
#
# Called from a `FindTPL<tplName>.cmake`_ module which first calls
# ``find_package(<externalPkg>)``and the calls this function to get and
# external package that uses modern CMake IMPORTED targets. This function
# creates the ``<tplName>::all_libs`` target and creates a TriBITS-compliant
# external package wrapper file `<tplName>Config.cmake`.
#
# Usage::
#
#   tribits_extpkg_create_imported_all_libs_target_and_config_file(
#     <tplName>
#     INNER_FIND_PACKAGE_NAME <externalPkg>
#     IMPORTED_TARGETS_FOR_ALL_LIBS <importedTarget0> <importedTarget1> ... )
#
# This function is called from a TriBITS ``FindTPL<tplName>.cmake`` wrapper
# module after it calls ``find_package(<externalPkg>)`` and then this function
# creates the IMPORTED target ``<tplName>::all_libs`` from the list of
# IMPORTED targets ``<importedTarget0> <importedTarget1> ...`` which are
# defined from the call ``find_package(<externalPkg>)``.  This function also
# takes care of generating the correct ``<tplName>Config.cmake`` file under
# the directory::
#
#   ${${PROJECT_NAME}_BINARY_DIR}/${${PROJECT_NAME}_BUILD_DIR_EXTERNAL_PKGS_DIR}
#
# The generated ``<tplName>Config.cmake`` file pulls in the upstream
# TriBITS-compliant ``<UpstreamPkg>Config.cmake` files, calls
# ``find_dependency(<externalPkg>)`` (with no other arguments), defines the
# `<tplName>::all_libs`` target, and then sets up the correct dependencies
# between these targets.
#
# For more details, see `Creating FindTPL<tplName>.cmake using find_package()
# with IMPORTED targets`_.
#
function(tribits_extpkg_create_imported_all_libs_target_and_config_file
    tplName
  )

  # Parse arguments
  cmake_parse_arguments(
     PARSE_ARGV 1
     PARSE "" "" # prefix, options, one_value_keywords
     "INNER_FIND_PACKAGE_NAME;IMPORTED_TARGETS_FOR_ALL_LIBS"  #multi_value_keywords
     )
  tribits_check_for_unparsed_arguments(PARSE)
  tribits_assert_parse_arg_one_value(PARSE  INNER_FIND_PACKAGE_NAME)
  tribits_assert_parse_arg_one_or_more_values(PARSE IMPORTED_TARGETS_FOR_ALL_LIBS)

  # Create imported target <tplName>::all_libs
  add_library(${tplName}::all_libs  INTERFACE  IMPORTED  GLOBAL)
  foreach (importedTarget  IN LISTS  PARSE_IMPORTED_TARGETS_FOR_ALL_LIBS)
    target_link_libraries(${tplName}::all_libs  INTERFACE  ${importedTarget})
  endforeach()

  # Create the TriBITS-compliant <tplName>Config.cmake wrapper file
  tribits_extpkgwit_create_package_config_file(
    ${tplName}
    INNER_FIND_PACKAGE_NAME ${PARSE_INNER_FIND_PACKAGE_NAME}
    IMPORTED_TARGETS_FOR_ALL_LIBS ${PARSE_IMPORTED_TARGETS_FOR_ALL_LIBS} )

endfunction()

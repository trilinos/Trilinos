# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER



# Macro that configures the package's main config.h file
#
function(tribits_add_config_define DEFINE)
  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("-- " "Package ${PARENT_PACKAGE_NAME}: adding compiler"
      " define to config file: ${DEFINE}")
  endif()
  global_set(${PARENT_PACKAGE_NAME}_CONFIG_DEFINES
    "${${PARENT_PACKAGE_NAME}_CONFIG_DEFINES}\n#define ${DEFINE}")
  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("-- ${${PARENT_PACKAGE_NAME}_CONFIG_DEFINES}")
  endif()
endfunction()


# @FUNCTION: tribits_configure_file()
#
# Macro that configures the package's main configured header file (typically
# called ``${PACKAGE_NAME}_config.h`` but any name can be used).
#
# Usage::
#
#   tribits_configure_file(<packageConfigFile>)
#
# This function requires the file::
#
#    ${PACKAGE_SOURCE_DIR}/cmake/<packageConfigFile>.in
#
# exists and it creates the file::
#
#   ${CMAKE_CURRENT_BINARY_DIR}/<packageConfigFile>
#
# by calling the built-in ``configure_file()`` command::
#
#   configure_file(
#     ${PACKAGE_SOURCE_DIR}/cmake/<packageConfigFile>.in
#     ${CMAKE_CURRENT_BINARY_DIR}/<packageConfigFile>
#     )
#
# which does basic substitution of CMake variables (see documentation for
# built-in CMake `configure_file()`_ command for rules on how it performs
# substitutions).  This command is typically used to configure the package's
# main `<packageDir>/cmake/<packageName>_config.h.in`_ file.
#
# In addition to just calling ``configure_file()``, this function also aids in
# creating configured header files adding macros for deprecating code as
# described below.
#
# **Deprecated Code Macros**
#
# If ``${PARENT_PACKAGE_NAME}_SHOW_DEPRECATED_WARNINGS`` is ``TRUE`` (see
# `tribits_add_show_deprecated_warnings_option()`_), then the local CMake
# variable ``${PARENT_PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS`` is set which
# adds a define ``<PARENT_PACKAGE_NAME_UC>_DEPRECATED`` (where
# ``<PARENT_PACKAGE_NAME_UC>`` is the package name in all upper-case letters)
# which adds a compiler-specific deprecated warning for an entity.  To take
# advantage of this, just add the line::
#
#   @<PARENT_PACKAGE_NAME_UC>_DEPRECATED_DECLARATIONS@
#
# to the ``<packageConfigFile>.in`` file and it will be expanded at configure
# time.
#
# Then C/C++ code can use this macro to deprecate functions, variables,
# classes, etc., for example, using::
#
#   <PARENT_PACKAGE_NAME_UC>_DEPRECATED class SomeDepreatedClass { ... }.
#
# If the particular compiler does not support deprecated warnings, then this
# macro is defined to be empty.  See `Regulated Backward Compatibility and
# Deprecated Code`_ for more details.
#
function(tribits_configure_file  PACKAGE_NAME_CONFIG_FILE)

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("\nPACKAGE_CONFIGURE_FILE: ${PACKAGE_NAME_CONFIG_FILE}")
  endif()

  # Set up the deprecated attribute if showing deprecated warnings
  if (${PARENT_PACKAGE_NAME}_SHOW_DEPRECATED_WARNINGS)
    multiline_set(${PARENT_PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS
      "#ifndef ${PARENT_PACKAGE_NAME_UC}_DEPRECATED\n"
      "#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))\n"
      "#    define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED  __attribute__((__deprecated__))\n"
      "#  else\n"
      "#    define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED\n"
      "#  endif\n"
      "#endif\n"
      "\n"
      "#ifndef ${PARENT_PACKAGE_NAME_UC}_DEPRECATED_MSG\n"
      "#  if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5))\n"
      "#    define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__ (#MSG) ))\n"
      "#  elif (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))\n"
      "#    define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__))\n"
      "#  else\n"
      "#    define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED_MSG(MSG)\n"
      "#  endif\n"
      "#endif\n"
      )
  else()
    multiline_set(${PARENT_PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS
      "#define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED\n"
      "#define ${PARENT_PACKAGE_NAME_UC}_DEPRECATED_MSG(MSG)\n"
      )
  endif()

  if (${PARENT_PACKAGE_NAME}_HIDE_DEPRECATED_CODE)
    string(APPEND ${PARENT_PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS
      "\n#define ${PARENT_PACKAGE_NAME_UC}_HIDE_DEPRECATED_CODE")
  endif()

  # Set up the macro to create the define for time monitor
  set(TIME_MONITOR_DEFINE_NAME ${PARENT_PACKAGE_NAME_UC}_TEUCHOS_TIME_MONITOR)
  set(FUNC_TIME_MONITOR_MACRO_NAME ${PARENT_PACKAGE_NAME_UC}_FUNC_TIME_MONITOR)
  set(FUNC_TIME_MONITOR_DIFF_MACRO_NAME ${PARENT_PACKAGE_NAME_UC}_FUNC_TIME_MONITOR_DIFF)
  if (${PARENT_PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR)
    multiline_set(${PARENT_PACKAGE_NAME_UC}_TEUCHOS_TIME_MONITOR_DECLARATIONS
      "#ifndef ${FUNC_TIME_MONITOR_MACRO_NAME}\n"
      "#  define ${TIME_MONITOR_DEFINE_NAME}\n"
      "#  define ${FUNC_TIME_MONITOR_MACRO_NAME}(FUNCNAME) \\\n"
      "     TEUCHOS_FUNC_TIME_MONITOR_DIFF(FUNCNAME, ${PARENT_PACKAGE_NAME_UC})\n"
      "#  define ${FUNC_TIME_MONITOR_DIFF_MACRO_NAME}(FUNCNAME, DIFF) \\\n"
      "     TEUCHOS_FUNC_TIME_MONITOR_DIFF(FUNCNAME, DIFF)\n"
      "#endif\n"
      )
  else()
    multiline_set(${PARENT_PACKAGE_NAME_UC}_TEUCHOS_TIME_MONITOR_DECLARATIONS
      "#define ${FUNC_TIME_MONITOR_MACRO_NAME}(FUNCNAME)\n"
      "#define ${FUNC_TIME_MONITOR_DIFF_MACRO_NAME}(FUNCNAME, DIFF)\n"
      )
  endif()

  # Configure the file
  configure_file(
    ${PACKAGE_SOURCE_DIR}/cmake/${PACKAGE_NAME_CONFIG_FILE}.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME_CONFIG_FILE}
    )

endfunction()

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


#
# This file gets included in the main TriBITS framework.  It is put here to
# reduce the size of the tribits/core/ directory.
#


function(tribits_write_deps_to_xml_string  PACKAGE_NAME  LIST_TYPE
  XML_VAR
  )

  set(LOC_XML "${${XML_VAR}}")

  set(DEPS_VAR ${PACKAGE_NAME}_${LIST_TYPE})
  assert_defined(DEPS_VAR)
  set(DEPS ${${DEPS_VAR}})

  #print_var(PACKAGE_NAME)
  #print_var(DEPS)

  if (NOT DEPS)

    list(APPEND LOC_XML
      "    <${LIST_TYPE}/>\n" )

  else()

    set(VALUE_STR "")

    foreach(DEP ${DEPS})

      if(VALUE_STR)
        set(VALUE_STR "${VALUE_STR},")
      endif()

      set(VALUE_STR "${VALUE_STR}${DEP}")

    endforeach()

    list(APPEND LOC_XML
      "    <${LIST_TYPE} value=\"${VALUE_STR}\"/>\n" )

  endif()

  if (LOC_XML)
    set(${XML_VAR} "${LOC_XML}" PARENT_SCOPE)
  endif()

endfunction()


#
# Function that writes the dependency information for ${PROJECT_NAME} into
# an XML file for other tools to use.
#

function(tribits_dump_deps_xml_file)

  set(DEPS_XM)

  get_filename_component(PROJECT_BASE_DIR_NAME  "${${PROJECT_NAME}_SOURCE_DIR}"  NAME)

  list(APPEND DEPS_XML
    "<PackageDependencies project=\"${PROJECT_NAME}\" baseDirName=\"${PROJECT_BASE_DIR_NAME}\">\n"
     )

  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES})

    #message("")
    #print_var(TRIBITS_PACKAGE)

    list(APPEND DEPS_XML
      "  <Package name=\"${TRIBITS_PACKAGE}\" dir=\"${${TRIBITS_PACKAGE}_REL_SOURCE_DIR}\" type=\"${${TRIBITS_PACKAGE}_TESTGROUP}\">\n")

    tribits_write_deps_to_xml_string(${TRIBITS_PACKAGE} LIB_REQUIRED_DEP_PACKAGES DEPS_XML)
    tribits_write_deps_to_xml_string(${TRIBITS_PACKAGE} LIB_OPTIONAL_DEP_PACKAGES DEPS_XML)
    tribits_write_deps_to_xml_string(${TRIBITS_PACKAGE} TEST_REQUIRED_DEP_PACKAGES DEPS_XML)
    tribits_write_deps_to_xml_string(${TRIBITS_PACKAGE} TEST_OPTIONAL_DEP_PACKAGES DEPS_XML)
    tribits_write_deps_to_xml_string(${TRIBITS_PACKAGE} LIB_REQUIRED_DEP_TPLS DEPS_XML)
    tribits_write_deps_to_xml_string(${TRIBITS_PACKAGE} LIB_OPTIONAL_DEP_TPLS DEPS_XML)
    tribits_write_deps_to_xml_string(${TRIBITS_PACKAGE} TEST_REQUIRED_DEP_TPLS DEPS_XML)
    tribits_write_deps_to_xml_string(${TRIBITS_PACKAGE} TEST_OPTIONAL_DEP_TPLS DEPS_XML)

    list(APPEND DEPS_XML
      "    <EmailAddresses>\n"
      "      <Regression address=\"${${TRIBITS_PACKAGE}_REGRESSION_EMAIL_LIST}\"/>\n"
      "    </EmailAddresses>\n"
      )

    list(APPEND DEPS_XML
      "    <ParentPackage value=\"${${TRIBITS_PACKAGE}_PARENT_PACKAGE}\"/>\n"
      )

    list(APPEND DEPS_XML
      "  </Package>\n" )

  endforeach()

  list(APPEND DEPS_XML
    "</PackageDependencies>\n" )

  #print_var(DEPS_XML)
  string(REPLACE "\n;" "\n" DEPS_XML "${DEPS_XML}")
  file(WRITE ${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE} "${DEPS_XML}" )

endfunction()


# @MACRO: tribits_write_xml_dependency_files()
#
# Usage::
#
#   tribits_write_xml_dependency_files()
#
# Macro that output XML dependency files if asked based in the global project
# package dependency graph previously constructed..
#
macro(tribits_write_xml_dependency_files)

  tribits_config_code_start_timer(WRITE_DEPENDENCY_FILES_TIME_START_SECONDS)

  #print_var(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE)
  if (${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE)
    if (NOT IS_ABSOLUTE ${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE})
      set(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE
        ${CMAKE_CURRENT_BINARY_DIR}/${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE})
    endif()
    message("" )
    message("Dumping the XML dependencies file ${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE} ..." )
    tribits_dump_deps_xml_file()
  endif()

  #print_var(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE)
  if (${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE AND ${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE)
    if (NOT IS_ABSOLUTE ${${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE})
      set(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE
        ${CMAKE_CURRENT_BINARY_DIR}/${${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE})
    endif()
    message("" )
    message("Dumping the HTML dependencies webpage file ${${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE} ..." )
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE}
        ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CI_SUPPORT_DIR}/dump-package-dep-table.py
        --input-xml-deps-file=${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE}
        --output-html-deps-file=${${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE} )
  endif()

  #print_var(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE)
  if (${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE AND ${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE)
    if (NOT IS_ABSOLUTE ${${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE})
      set(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/${${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE})
    endif()
    message("" )
    message("Dumping the CDash XML dependencies file ${${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE} ..." )
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE}
        ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CTEST_DRIVER_DIR}/dump-cdash-deps-xml-file.py
        --input-xml-deps-file=${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE}
        --output-cdash-deps-xml-file=${${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE} )
  endif()

  if (${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE
    OR ${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE
    OR ${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE
    )
    tribits_config_code_stop_timer(WRITE_DEPENDENCY_FILES_TIME_START_SECONDS
      "\nTotal time to write dependency files")
  endif()

endmacro()

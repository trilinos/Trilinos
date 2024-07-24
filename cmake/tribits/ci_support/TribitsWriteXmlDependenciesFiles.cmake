# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


################################################################################
# NOTE: This file gets included in the main TriBITS core framework.  It is put
# here to reduce the size of the tribits/core/ directory as this is optional
# behavior.
################################################################################


include(TribitsConfigureTiming)


# @MACRO: tribits_write_xml_dependency_files()
#
# Usage::
#
#   tribits_write_xml_dependency_files()
#
# Macro that outputs XML dependency files if asked based in the global project
# package dependency graph previously constructed.
#
macro(tribits_write_xml_dependency_files)

  tribits_config_code_start_timer(WRITE_DEPENDENCY_FILES_TIME_START_SECONDS)

  if (${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE)
    if (NOT IS_ABSOLUTE ${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE})
      set(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE
        ${CMAKE_CURRENT_BINARY_DIR}/${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE})
    endif()
    message("" )
    message("Dumping the XML dependencies file"
      " ${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE} ..." )
    tribits_dump_deps_xml_file()
  endif()

  if (${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE AND ${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE)
    if (NOT IS_ABSOLUTE ${${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE})
      set(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE
        ${CMAKE_CURRENT_BINARY_DIR}/${${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE})
    endif()
    message("" )
    message("Dumping the HTML dependencies webpage file"
      " ${${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE} ..." )
    set(tribitsCiSupportDir "${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CI_SUPPORT_DIR}")
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE}
        ${tribitsCiSupportDir}/dump-package-dep-table.py
        --input-xml-deps-file=${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE}
        --output-html-deps-file=${${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE} )
  endif()

  if (${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE AND ${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE)
    if (NOT IS_ABSOLUTE ${${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE})
      set(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE
	${CMAKE_CURRENT_BINARY_DIR}/${${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE})
    endif()
    message("" )
    message("Dumping the CDash XML dependencies file"
      " ${${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE} ..." )
    set(tribitsCtestDriverDir
      "${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CTEST_DRIVER_DIR}")
    if (EXISTS "${tribitsCtestDriverDir}")
      execute_process(
        COMMAND ${PYTHON_EXECUTABLE}
          ${tribitsCtestDriverDir}/dump-cdash-deps-xml-file.py
          --input-xml-deps-file=${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE}
          --output-cdash-deps-xml-file=${${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE})
    else()
      message(FATAL_ERROR "\nERROR: Can't write"
	" ${${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE}"
	" because '${tribitsCtestDriverDir}' does not exist!")
    endif()
  endif()

  if (${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE
      OR ${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE
      OR ${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE
    )
    tribits_config_code_stop_timer(WRITE_DEPENDENCY_FILES_TIME_START_SECONDS
      "\nTotal time to write dependency files")
  endif()

endmacro()


# Function that writes the dependency information for the TriBITS project
# ${PROJECT_NAME} into an XML file for other tools to use.
#
function(tribits_dump_deps_xml_file)

  set(depsXml "")

  get_filename_component(PROJECT_BASE_DIR_NAME  "${${PROJECT_NAME}_SOURCE_DIR}"  NAME)

  list(APPEND depsXml
    "<PackageDependencies project=\"${PROJECT_NAME}\" baseDirName=\"${PROJECT_BASE_DIR_NAME}\">\n"
     )

  foreach(tribitsPackage ${${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES})

    list(APPEND depsXml
      "  <Package name=\"${tribitsPackage}\" dir=\"${${tribitsPackage}_REL_SOURCE_DIR}\" type=\"${${tribitsPackage}_TESTGROUP}\">\n")

    tribits_write_deps_to_xml_string(${tribitsPackage} LIB REQUIRED PACKAGES depsXml)
    tribits_write_deps_to_xml_string(${tribitsPackage} LIB OPTIONAL PACKAGES depsXml)
    tribits_write_deps_to_xml_string(${tribitsPackage} TEST REQUIRED PACKAGES depsXml)
    tribits_write_deps_to_xml_string(${tribitsPackage} TEST OPTIONAL PACKAGES depsXml)
    tribits_write_deps_to_xml_string(${tribitsPackage} LIB REQUIRED TPLS depsXml)
    tribits_write_deps_to_xml_string(${tribitsPackage} LIB OPTIONAL TPLS depsXml)
    tribits_write_deps_to_xml_string(${tribitsPackage} TEST REQUIRED TPLS depsXml)
    tribits_write_deps_to_xml_string(${tribitsPackage} TEST OPTIONAL TPLS depsXml)

    list(APPEND depsXml
      "    <EmailAddresses>\n"
      "      <Regression address=\"${${tribitsPackage}_REGRESSION_EMAIL_LIST}\"/>\n"
      "    </EmailAddresses>\n"
      )

    list(APPEND depsXml
      "    <ParentPackage value=\"${${tribitsPackage}_PARENT_PACKAGE}\"/>\n"
      )

    list(APPEND depsXml
      "  </Package>\n" )

  endforeach()

  list(APPEND depsXml
    "</PackageDependencies>\n" )

  string(REPLACE "\n;" "\n" depsXml "${depsXml}")
  file(WRITE ${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE} "${depsXml}" )

endfunction()


function(tribits_write_deps_to_xml_string  packageName libOrTest  requiredOrOptional
    packagesOrTpls  xmlVarInOut
  )

  set(localXml "${${xmlVarInOut}}")

  set(listType ${libOrTest}_${requiredOrOptional}_DEP_${packagesOrTpls})
  message("")

  tribits_get_legacy_package_deps_sublist(${packageName} ${libOrTest}
    ${requiredOrOptional} ${packagesOrTpls} legacyPackageDepsList)

  if (NOT legacyPackageDepsList)
    list(APPEND localXml
      "    <${listType}/>\n" )
  else()
    set(depsListStr "")
    foreach(depPkg IN LISTS legacyPackageDepsList)
      tribits_append_dep_to_xml_string(${depPkg} depsListStr)
    endforeach()
    list(APPEND localXml
      "    <${listType} value=\"${depsListStr}\"/>\n" )
  endif()

  if (localXml)
    set(${xmlVarInOut} "${localXml}" PARENT_SCOPE)
  endif()

endfunction()


function(tribits_get_legacy_package_deps_sublist  packageName  libOrTest
    requiredOrOptional  packagesOrTpls  legacyPackageDepsListOut
  )

  set(legacyPackageDepsList "")

  foreach(depPkg IN LISTS ${packageName}_${libOrTest}_DEFINED_DEPENDENCIES)

    set(matchesRequriedOrOptional FALSE)
    if (((requiredOrOptional STREQUAL "REQUIRED")
        AND ${packageName}_${libOrTest}_DEP_REQUIRED_${depPkg})
      OR
        ((requiredOrOptional STREQUAL "OPTIONAL")
        AND (NOT ${packageName}_${libOrTest}_DEP_REQUIRED_${depPkg}))
      )
      set(matchesRequriedOrOptional TRUE)
    endif()

    set(matchesPackagesOrTpls FALSE)
    if (((packagesOrTpls STREQUAL "PACKAGES")
        AND (${depPkg}_PACKAGE_BUILD_STATUS STREQUAL "INTERNAL"))
      OR
        ((packagesOrTpls STREQUAL "TPLS")
        AND (${depPkg}_PACKAGE_BUILD_STATUS STREQUAL "EXTERNAL"))
      )
      set(matchesPackagesOrTpls TRUE)
    endif()

    if (matchesRequriedOrOptional AND matchesPackagesOrTpls)
      list(APPEND legacyPackageDepsList ${depPkg})
    endif()

  endforeach()

  set(${legacyPackageDepsListOut} "${legacyPackageDepsList}" PARENT_SCOPE)

endfunction()



function(tribits_append_dep_to_xml_string  depPkg  depsListStrInOut)
  set(depsListStr "${${depsListStrInOut}}")
  if (depsListStr)
    string(APPEND depsListStr ",")
  endif()
  string(APPEND depsListStr "${depPkg}")
  set(${depsListStrInOut} "${depsListStr}" PARENT_SCOPE)
endfunction() 

# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


#
# Top-level project logic to generate resources spec file
#
function(tribits_generate_ctest_resource_spec_file_project_logic)
  if (${PROJECT_NAME}_AUTOGENERATE_TEST_RESOURCE_FILE)
    if (CTEST_RESOURCE_SPEC_FILE STREQUAL CTEST_RESOURCE_SPEC_FILE_DEFAULT)
      tribits_generate_ctest_resource_spec_file()
    else()
      message("NOTE: The test resource file CTEST_RESOURCE_SPEC_FILE='${CTEST_RESOURCE_SPEC_FILE}'"
        " will not be auto-generated even through"
        " ${PROJECT_NAME}_AUTOGENERATE_TEST_RESOURCE_FILE=${${PROJECT_NAME}_AUTOGENERATE_TEST_RESOURCE_FILE}"
        " because its location does not match the default"
        " location '${CTEST_RESOURCE_SPEC_FILE_DEFAULT}'."
        "  If you want to auto-generate this file, please clear CTEST_RESOURCE_SPEC_FILE and"
        " reconfigure or create that file on your own and clear"
        " ${PROJECT_NAME}_AUTOGENERATE_TEST_RESOURCE_FILE."
        )
    endif()
  endif()
endfunction()


#
# Generate resource spec file
#
function(tribits_generate_ctest_resource_spec_file)
  set(GPUS_JSON)
  math(EXPR LAST_GPU "${${PROJECT_NAME}_CUDA_NUM_GPUS} - 1")
  set(FIRST 1)
  foreach(GPU RANGE 0 ${LAST_GPU})
    if(NOT FIRST)
      string(APPEND GPUS_JSON ",\n")
    endif()
    set(FIRST 0)
    string(APPEND GPUS_JSON "        {
          \"id\": \"${GPU}\",
          \"slots\": ${${PROJECT_NAME}_CUDA_SLOTS_PER_GPU}
        }")
  endforeach()
  file(WRITE "${CMAKE_BINARY_DIR}/ctest_resources.json" "{
  \"version\": {
    \"major\": 1,
    \"minor\": 0
  },
  \"local\": [
    {
      \"gpus\": [
${GPUS_JSON}
      ]
    }
  ]
}
")
endfunction()

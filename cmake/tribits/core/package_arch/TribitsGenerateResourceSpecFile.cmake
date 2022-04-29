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

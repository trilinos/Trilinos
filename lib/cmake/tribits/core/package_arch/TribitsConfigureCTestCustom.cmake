# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(TribitsGeneralMacros)

function(tribits_configure_ctest_custom  TRIBITS_PROJECT_SOURCE_DIR  OUTPUT_BINARY_DIR)
  set(CTEST_CUSTOM_IN ${TRIBITS_PROJECT_SOURCE_DIR}/cmake/ctest/CTestCustom.cmake.in)
  #print_var(CTEST_CUSTOM_IN)
  if(EXISTS ${CTEST_CUSTOM_IN})
    tribits_trace_file_processing(PROJECT  CONFIGURE "${CTEST_CUSTOM_IN}")
    configure_file(${CTEST_CUSTOM_IN} ${OUTPUT_BINARY_DIR}/CTestCustom.cmake)
  endif()
endfunction()

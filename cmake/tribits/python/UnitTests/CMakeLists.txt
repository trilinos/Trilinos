# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
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



INCLUDE(TribitsAddAdvancedTest)

# The mock test project 
SET(MOCK_PROJECT_NAME Trilinos)


SET(MOCK_TRILINOS_DIR ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_MOCK_TRILINOS_DIR})


IF (WIN23 AND NOT CYGWIN)
  SET(DIFF_COMMAND_STUFF "${CMAKE_COMMAND}" ARGS -E compare_files)
ELSE()
  SET(DIFF_COMMAND_STUFF diff ARGS)
ENDIF()


TRIBITS_ADD_ADVANCED_TEST( XmlHtmlCDashOutputTest
  OVERALL_WORKING_DIRECTORY TEST_NAME
  TEST_0 CMND ${CMAKE_COMMAND}
    ARGS
      -D${MOCK_PROJECT_NAME}_TRIBITS_DIR:FILEPATH=${${PROJECT_NAME}_TRIBITS_DIR}
      -D ${MOCK_PROJECT_NAME}_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING:BOOL=ON
      -D ${MOCK_PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=ON
      -D ${MOCK_PROJECT_NAME}_DEPS_XML_OUTPUT_FILE:STRING=TrilinosPackageDependencies.xml
      -D ${MOCK_PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE:STRING=TrilinosPackageDependenciesTable.html
      -D ${MOCK_PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE:STRING=CDashSubprojectDependencies.xml
      -D ${MOCK_PROJECT_NAME}_ENABLE_C:BOOL=OFF
      -D ${MOCK_PROJECT_NAME}_ENABLE_CXX:BOOL=OFF
      -D ${MOCK_PROJECT_NAME}_ENABLE_Fortran:BOOL=OFF
      ${MOCK_TRILINOS_DIR}
    PASS_REGULAR_EXPRESSION "Shortcircuiting after dependency tracking ..."
  TEST_1 CMND ${DIFF_COMMAND_STUFF}
    ${CMAKE_CURRENT_SOURCE_DIR}/TrilinosPackageDependencies.gold.xml
    TrilinosPackageDependencies.xml
  TEST_2 CMND ${DIFF_COMMAND_STUFF}
    ${CMAKE_CURRENT_SOURCE_DIR}/TrilinosPackageDependenciesTable.gold.html
     TrilinosPackageDependenciesTable.html
  TEST_3 CMND ${DIFF_COMMAND_STUFF}
    ${CMAKE_CURRENT_SOURCE_DIR}/CDashSubprojectDependencies.gold.xml
    CDashSubprojectDependencies.xml
  )

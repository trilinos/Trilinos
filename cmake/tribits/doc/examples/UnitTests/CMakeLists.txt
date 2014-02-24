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


TRIBITS_ADD_ADVANCED_TEST( RawHelloWorld
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  TEST_0 CMND ${CMAKE_COMMAND}
    ARGS
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/RawHelloWorld
    PASS_REGULAR_EXPRESSION_ALL
      "Configuring done"
      "Generating done"
      "Build files have been written to: .*tribits/doc/examples/UnitTests/TriBITS_RawHelloWorld"
  TEST_1 CMND make
    ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Built target hello_world_lib"
      "Built target hello_world"
      "Built target unit_tests"
  TEST_2 CMND ${CMAKE_CTEST_COMMAND}
    PASS_REGULAR_EXPRESSION_ALL
      ": test .*  Passed"
      ": unit_test .*  Passed"
      "100% tests passed, 0 tests failed out of 2"
  CATEGORIES NIGHTLY
  )


TRIBITS_ADD_ADVANCED_TEST( TribitsHelloWorld
  OVERALL_WORKING_DIRECTORY TEST_NAME
  OVERALL_NUM_MPI_PROCS 1
  TEST_0 CMND ${CMAKE_COMMAND}
    ARGS
      ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsHelloWorld
    PASS_REGULAR_EXPRESSION_ALL
      "Configuring done"
      "Generating done"
      "Build files have been written to: .*tribits/doc/examples/UnitTests/TriBITS_TribitsHelloWorld"
  TEST_1 CMND make
    ARGS ${CTEST_BUILD_FLAGS}
    PASS_REGULAR_EXPRESSION_ALL
      "Built target hello_world_lib"
      "Built target hello_world"
      "Built target HelloWorld_unit_tests"
  TEST_2 CMND ${CMAKE_CTEST_COMMAND}
    PASS_REGULAR_EXPRESSION_ALL
      ": HelloWorld_hello_world .*  Passed"
      ": HelloWorld_unit_tests .*  Passed"
      "100% tests passed, 0 tests failed out of 2"
  CATEGORIES NIGHTLY
  )

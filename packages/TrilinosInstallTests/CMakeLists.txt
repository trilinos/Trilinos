# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
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
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER

tribits_package(TrilinosInstallTests)


tribits_add_advanced_test(doInstall
  OVERALL_NUM_MPI_PROCS 1

  TEST_0
    MESSAGE "Removing installation directory if it already exists"
    CMND "${CMAKE_COMMAND}"
    ARGS "-DDIR_TO_REMOVE=${PROJECT_BINARY_DIR}/install"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/remove_dir_if_exists.cmake"

  TEST_1
    MESSAGE "Install whatever Trilinos packages have been enabled"
    CMND "${CMAKE_COMMAND}"
    ARGS --install ${PROJECT_BINARY_DIR}
      --prefix ${PROJECT_BINARY_DIR}/install
    OUTPUT_FILE  doInstall.out
    NO_ECHO_OUTPUT

  ADDED_TEST_NAME_OUT  doInstall_name
  )
  # NOTE: Above works even if Trilinos was configured without setting
  # -DCMAKE_PREFIX_PATH=<prefix> and tests proper usage of the install()
  # command.  However, note that it is not a perfect installation test because
  # the source dir and the build dir will still be sticking around in the
  # below example build.


tribits_add_advanced_test(find_package_Trilinos
  OVERALL_NUM_MPI_PROCS 1
  TEST_0
    MESSAGE "Configure a dummy project that calls find_package(Trilinos)"
    CMND "${CMAKE_COMMAND}"
    ARGS
      -DCMAKE_PREFIX_PATH=${PROJECT_BINARY_DIR}/install
      ${CMAKE_CURRENT_SOURCE_DIR}/find_package_Trilinos
    PASS_REGULAR_EXPRESSION_ALL
      "Trilinos_FOUND = '1'"
      "-- Configuring done"
      "-- Generating done"
    ALWAYS_FAIL_ON_NONZERO_RETURN
  ADDED_TEST_NAME_OUT  find_package_Trilinos_name
  )
  # NOTE: The above test will run find_package(Trilinos) for **all** of the
  # enabled packages, not just those needed by the project
  # simpleBuildAgainstTrilinos below.  This will show if calling
  # find_package(Trilinos) has any problems itself.

if (find_package_Trilinos_name)
  set_tests_properties(${find_package_Trilinos_name}
    PROPERTIES DEPENDS ${doInstall_name} )
endif()


tribits_add_advanced_test(simpleBuildAgainstTrilinos
  OVERALL_NUM_MPI_PROCS 1
  OVERALL_WORKING_DIRECTORY TEST_NAME
  EXCLUDE_IF_NOT_TRUE  ${PROJECT_NAME}_ENABLE_Tpetra  TPL_ENABLE_MPI

  TEST_0
    MESSAGE "Configure simpleBuildAgainstTrilinos against local Trilinos installation"
    CMND "${CMAKE_COMMAND}"
    ARGS
      -DCMAKE_PREFIX_PATH=${PROJECT_BINARY_DIR}/install
      -DCMAKE_BUILD_TYPE=DEBUG
      ${PROJECT_SOURCE_DIR}/demos/simpleBuildAgainstTrilinos
    PASS_REGULAR_EXPRESSION_ALL
      "-- Configuring done"
      "-- Generating done"
    ALWAYS_FAIL_ON_NONZERO_RETURN

  TEST_1
    MESSAGE "Build simpleBuildAgainstTrilinos"
    CMND "${CMAKE_COMMAND}"
    ARGS --build .

  TEST_2
    MESSAGE "Test simpleBuildAgainstTrilinos"
    CMND "${CMAKE_CTEST_COMMAND}"
    ARGS -VV
    PASS_REGULAR_EXPRESSION_ALL
      "Test #1: MyAppTest [.]+ +Passed"
      "100% tests passed, 0 tests failed out of 1"
    ALWAYS_FAIL_ON_NONZERO_RETURN

  ADDED_TEST_NAME_OUT  simpleBuildAgainstTrilinos_name
  )

if (simpleBuildAgainstTrilinos_name)
  set_tests_properties(${simpleBuildAgainstTrilinos_name}
    PROPERTIES DEPENDS ${doInstall_name} )
endif()


tribits_package_postprocess()


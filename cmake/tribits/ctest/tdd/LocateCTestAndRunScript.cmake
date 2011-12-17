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

# This script is a cmake -P script that should be called with the following
# variables defined:
#
#    -D binary_dir=${CMAKE_CURRENT_BINARY_DIR}
#    -D source_dir=${CMAKE_CURRENT_SOURCE_DIR}
#    -D ctest_type=${ctest_type}
#    -D scriptname=${scriptname}
#    -D TD_BASE_DIR=${TD_BASE_DIR}
#    -D testname=${testname}
#
# It looks recursively under ${TD_BASE_DIR}/tools/cmake-${ctest_type}
# for a ctest executable and uses it to drive a ctest -S script to run
# a dashboard for the project. It has to be run this way indirectly
# through a cmake -P script because the desired ctest executable
# location is unknown at CMake configure time of the driver project.

if(NOT CTEST_EXE)
  if(WIN32)
    set(ctest_filename "ctest.exe")
  else()
    set(ctest_filename "ctest")
  endif()

  message("globbing for '${ctest_filename}' of type '${ctest_type}'...")

  file(GLOB_RECURSE CTEST_EXE
    "${TD_BASE_DIR}/tools/cmake-${ctest_type}/${ctest_filename}")
endif()

if(NOT CTEST_EXE)
  message(FATAL_ERROR "error: '${ctest_type}' ctest could not be found...")
endif()

if(NOT EXISTS "${CTEST_EXE}")
  message(FATAL_ERROR "error: CTEST_EXE='${CTEST_EXE}' does not exist...")
endif()

message("CTEST_EXE='${CTEST_EXE}'")
execute_process(COMMAND ${CTEST_EXE} --version)

message("=========== variables ===========")
message("binary_dir='${binary_dir}'")
message("source_dir='${source_dir}'")
message("ctest_type='${ctest_type}'")
message("scriptname='${scriptname}'")
message("TD_BASE_DIR='${TD_BASE_DIR}'")
message("testname='${testname}'")
message("=================================")

message("========== environment ==========")
execute_process(COMMAND ${CMAKE_COMMAND} -E environment)
message("=================================")

message("============ script =============")
message("executing ctest -S '${scriptname}' for test '${testname}'...")

execute_process(COMMAND ${CTEST_EXE}
  -S
  "${source_dir}/${scriptname}"
  -V
  --output-log
  "${binary_dir}/${testname}.log"
  RESULT_VARIABLE rv
  )

if(NOT "${rv}" STREQUAL "0")
  message("warning: calling ctest -S script failed with '${rv}'")
endif()

message("=================================")

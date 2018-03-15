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


INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../../TrilinosCTestDriverCore.cmake")

# Must set the Jenkins JOB_NAME which will be the CDash build name 
IF ("$ENV{JOB_NAME}" STREQUAL "")
  MESSAGE(FATAL_ERROR "Error, must set env var JOB_NAME")
ENDIF()
SET(CTEST_BUILD_NAME "$ENV{JOB_NAME}")
SET(CTEST_TEST_TIMEOUT 600)
SET( CTEST_NOTES_FILES
  "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}"
  "${TRIBITS_PROJECT_ROOT}/cmake/std/sems/atdm/load_atdm_7.2_dev_env.sh"
  )
SET( CTEST_BUILD_FLAGS "-j10 -i" )
SET( CTEST_PARALLEL_LEVEL "10" )
SET(Trilinos_ENABLE_SECONDARY_TESTED_CODE OFF)
SET(Trilinos_ENABLE_CONFIGURE_TIMING ON)
SET(Trilinos_BRANCH develop)
SET(Trilinos_EXTRAREPOS_FILE NONE)
# Just test all the packages
#SET(Trilinos_PACKAGES)
SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/sems/atdm/SEMSATDMSettings.cmake,cmake/std/MpiReleaseDebugSharedPtSettings.cmake,cmake/std/BasicCiTestingSettings.cmake,cmake/std/sems/SEMSDevEnv.cmake"
  "-DTrilinos_ENABLE_CONFIGURE_TIMING=ON"
  )
SET(CTEST_TEST_TYPE Nightly)
SET(Trilinos_TRACK ATDM)

# Must set the site name so that it does not change depending on what node
# runs the build.
SET(CTEST_SITE "sems-rhel6")

TRIBITS_CTEST_DRIVER()

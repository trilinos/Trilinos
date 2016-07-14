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

# TEMPORARY SETUP FOR TESTING

MESSAGE( "Turning CTEST_DO_UPDATES off to prevent attempt to access Sandia libraries." )
SET(CTEST_DO_UPDATES FALSE)

MESSAGE( "Turning CTEST_DO_SUBMIT off to prevent dashboard submissions temporarily." )
SET(CTEST_DO_SUBMIT FALSE)
SET(CTEST_SUBMIT_CDASH_SUBPROJECTS_DEPS_FILE FALSE)


#
# Set the options specific to this build case - must be BEFORE including the jenkins file
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE DEBUG)
#SET(BUILD_DIR_NAME)
SET(CTEST_PARALLEL_LEVEL 8)
SET(CTEST_TEST_TYPE Experimental)
SET(CTEST_TEST_TIMEOUT 900)
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)

# Generic driver script for building & testing Trilinos with Jenkins
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.generic.jenkins.cmake")

SET(Trilinos_PACKAGES Zoltan2)

# The Trilinos

SET(EXTRA_CONFIGURE_OPTIONS
"-DZoltan2_ENABLE_Experimental:BOOL=$ENV{JENKINS_Zoltan2_ENABLE_Experimental}"
"-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=$ENV{JENKINS_ENABLE_EXPLICIT_INSTANTIATION}"
"-DTeuchos_ENABLE_LONG_LONG_INT:BOOL=$ENV{JENKINS_Teuchos_ENABLE_LONG_LONG_INT}"

"-DTpetra_INST_INT_INT:BOOL=$ENV{JENKINS_Tpetra_INST_INT_INT}"
"-DTpetra_INST_INT_UNSIGNED:BOOL=$ENV{JENKINS_Tpetra_INST_INT_UNSIGNED}"
"-DTpetra_INST_INT_LONG:BOOL=$ENV{JENKINS_Tpetra_INST_INT_LONG}"
"-DTpetra_INST_INT_LONG_LONG:BOOL=$ENV{JENKINS_Tpetra_INST_INT_LONG_LONG}"
"-DTpetra_INST_INT_UNSIGNED_LONG:BOOL=$ENV{JENKINS_Tpetra_INST_INT_UNSIGNED_LONG}"

"-DTpetra_INST_DOUBLE:BOOL=$ENV{JENKINS_Tpetra_INST_DOUBLE}"
"-DTpetra_INST_FLOAT:BOOL=$ENV{JENKINS_Tpetra_INST_FLOAT}"
"-DTpetra_INST_COMPLEX_DOUBLE:BOOL=$ENV{JENKINS_Tpetra_INST_COMPLEX_DOUBLE}"
"-DTpetra_INST_COMPLEX_FLOAT:BOOL=$ENV{JENKINS_Tpetra_INST_COMPLEX_FLOAT}"

"-DTPL_ENABLE_SCOTCH:BOOL=$ENV{JENKINS_ENABLE_SCOTCH}"

)

# SET(EXTRA_EXCLUDE_PACKAGES MOOCHO Sundance CTrilinos ForTrilinos Optika Mesquite preCopyrightTrilinos TerminalApplication)

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

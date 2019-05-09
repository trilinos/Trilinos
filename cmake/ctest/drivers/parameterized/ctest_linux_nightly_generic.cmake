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


INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.generic.cmake")

#
# Set the options specific to this build case
#
SET(COMM_TYPE $ENV{JENKINS_COMM_TYPE})
SET(BUILD_TYPE $ENV{JENKINS_BUILD_TYPE})


STRING(REPLACE "sems-" "" COMPILER_DIR $ENV{COMPILER_MODULE})
STRING(REPLACE "/" "_"    COMPILER_DIR ${COMPILER_DIR})


SET(COMM_DIR "")
IF(COMM_TYPE STREQUAL MPI)
    STRING(REPLACE "sems-" "" COMM_DIR $ENV{MPI_MODULE})
    STRING(REPLACE "/" "_"    COMM_DIR ${COMM_DIR})
ENDIF()


# If JENKINS_BUILD_TAG is nonempty, prepend an underscore to it, otherwise leave
# it blank.  If JENKINS_BUILD_TAG doesn't exist, it shoud be set to empty.
SET(BUILD_TAG "$ENV{JENKINS_BUILD_TAG}")
IF(BUILD_TAG)
    # You'd think we could just use "^_*" as the pattern and "_" as the replace
    # but CMake isn't happy with that.
    STRING(REGEX REPLACE "^_+" "" BUILD_TAG ${BUILD_TAG})
    SET(BUILD_TAG "_${BUILD_TAG}")
ENDIF()


# Get the Jenkins Environment variable, if it exists, containing extra parameters 
# for EXTRA_CONFIGURE_OPTIONS
# The string should come in like looking like this:
#    '"-DTrilinos_ENABLE_FOO:BOOL=ON" "-DTrilinos_ENABLE_BAR:BOOL=ON"'
#
SET(JENKINS_EXTRA_CONFIGURE_OPTIONS "$ENV{JENKINS_Trilinos_EXTRA_CONFIGURE_OPTIONS}")
IF(JENKINS_EXTRA_CONFIGURE_OPTIONS)
    MESSAGE(STATUS "JENKINS_Trilinos_EXTRA_CONFIGURE_OPTIONS: ${JENKINS_EXTRA_CONFIGURE_OPTIONS}")

    # Strip leading and trailing whitespace
    STRING(STRIP "${JENKINS_EXTRA_CONFIGURE_OPTIONS}" JENKINS_EXTRA_CONFIGURE_OPTIONS)

    # There should just be one leading and trailing '"' character,
    # But when JENKINS_EXTRA_CONFIGURE_OPTIONS is added to EXTRA_CONFIGURE_OPTIONS 
    # CMake is adding in an extra quotation by the time this gets to the command line.
    # So let's strip off the leading and trailing quotation from the string itself.
    STRING(REGEX REPLACE "^\"+" "" JENKINS_EXTRA_CONFIGURE_OPTIONS ${JENKINS_EXTRA_CONFIGURE_OPTIONS})
    STRING(REGEX REPLACE "\"+$" "" JENKINS_EXTRA_CONFIGURE_OPTIONS ${JENKINS_EXTRA_CONFIGURE_OPTIONS})
ENDIF()

# Debugging Messages
# MESSAGE("")
# MESSAGE("---")
# MESSAGE("--- ENV:JENKINS_EXTRA_CONFIGURE_OPTIONS: '$ENV{JENKINS_Trilinos_EXTRA_CONFIGURE_OPTIONS}'")
# MESSAGE("--- JENKINS_EXTRA_CONFIGURE_OPTIONS....: '${JENKINS_EXTRA_CONFIGURE_OPTIONS}'")
# MESSAGE("--- BUILD_TAG..........................: '${BUILD_TAG}'")
# MESSAGE("---")

# Set the build name of the job (this is reported to CDash as the job name)
SET(BUILD_DIR_NAME ${COMM_TYPE}_${BUILD_TYPE}_${COMPILER_DIR}_${COMM_DIR}${BUILD_TAG}_DEV)

SET(CTEST_PARALLEL_LEVEL 16)


# Note: CTEST_TEST_TYPE drives some side-effects in Tribits that should be 
#       taken into account.  If CTEST_TEST_TYPE is Experimental, Tribits will
#       override Trilinos_TRACK and *always* assign results to the Experimental
#       track on CDash.  Also, Tribits does different things based on CTEST_TEST_TYPE
#       being one of the 3 blessed types (Nightly, Continuous, Experimental), so it's
#       best to keep CTEST_TEST_TYPE one of these values.  
SET(CTEST_TEST_TYPE Nightly)
SET(Trilinos_TRACK  $ENV{JENKINS_JOB_TYPE})


SET(CTEST_TEST_TIMEOUT 900)


SET(EXTRA_CONFIGURE_OPTIONS
    "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
    "-DTPL_ENABLE_SuperLU:BOOL=ON"
    "-DTPL_ENABLE_Netcdf:BOOL=ON"
    "-DTPL_ENABLE_HDF5:BOOL=ON"
    "-DTPL_ENABLE_Boost:BOOL=ON"
    "-DTPL_ENABLE_BoostLib:BOOL=ON"
    "-DTPL_ENABLE_Zlib:BOOL=ON"

    "-DMueLu_INST_DOUBLE_INT_LONGINT:BOOL=ON"
    "-DTPL_ENABLE_Matio:BOOL=OFF"

    # wcmclen: 2017-07-10: Disabling tests Anasazi_Epetra_OrthoManagerGenTester_1_MPI_4 and
    #                      Anasazi_Epetra_ModalSolversTester_MPI_4 due to unstable behaviour
    #                      in testing.  Heidi has been looking into it but it appears that
    #                      it's not an easy fix.
    #                      Reference: Trilinos issue #1393 (https://github.com/trilinos/Trilinos/issues/1393)
    "-DAnasazi_Epetra_ModalSolversTester_MPI_4_DISABLE:BOOL=ON"
    "-DAnasazi_Epetra_OrthoManagerGenTester_0_MPI_4_DISABLE:BOOL=ON"
    "-DAnasazi_Epetra_OrthoManagerGenTester_1_MPI_4_DISABLE:BOOL=ON"

    # Append extra options from the Jenkins parameter JENKINS_EXTRA_CONFIGURE_OPTIONS
    ${JENKINS_EXTRA_CONFIGURE_OPTIONS}
)

#"-DMPI_EXEC_POST_NUMPROCS_FLAGS:STRING=-bind-to;socket;--map-by;socket"
#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

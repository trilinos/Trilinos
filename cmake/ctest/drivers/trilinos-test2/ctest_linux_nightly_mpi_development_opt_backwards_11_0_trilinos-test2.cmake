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

#run command to create a tarball and install from that so we can test against that installation

#hard coding where the installation dir is for now.
#Ideally this should be pulled in from the same variable that is used to set the installation
#dir for the script.
SET(INSTALLATION_DIR "${CMAKE_CURRENT_BINARY_DIR}/../../installation/installed")

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.trilinos-test2.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME MPI_OPT_DEV_BACKWARDS_COMPATIBILITY_11.0)
#SET(CTEST_TEST_TIMEOUT 900)
SET(CTEST_TEST_TYPE EXPERIMENTAL)

#setting this temporarily while this is an experimental test
#the default repository for experimental is /space/git/Trilinos
#but the branch in question only exists on /space/git/nightly/Trilinos
SET(Trilinos_REPOSITORY_LOCATION "software.sandia.gov:/space/git/nightly/Trilinos")

SET(Trilinos_BRANCH "trilinos-release-11-0-branch")

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE OFF)

SET(EXTRA_EXCLUDE_PACKAGES Mesquite RBGen)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_INSTALLATION_DIR=${INSTALLATION_DIR}"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-DMPI_BASE_DIR:PATH=/home/trilinos"
  "-DTrilinos_ENABLE_DEVELOPMENT_MODE:BOOL=OFF"
  "-DTPL_ENABLE_Pthread:BOOL=ON"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  "-DBoost_INCLUDE_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/boost-1.49.0"
  "-DNetcdf_LIBRARY_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/pnetcdf_4.2/lib"
  "-DNetcdf_INCLUDE_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/pnetcdf_4.2/include"
  "-DHDF5_INCLUDE_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/phdf5-1.8.6/include"
  "-DHDF5_LIBRARY_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/phdf5-1.8.6/lib"
  "-DTrilinos_ENABLE_RBGen=OFF"
  "-DTPL_ENABLE_SuperLU=OFF"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

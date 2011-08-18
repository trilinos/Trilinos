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


INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.zan.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME MPI_OPT_DEV2)
SET(COMPILER_VERSION "GCC-4.6.1")
SET(ENV{LD_LIBRARY_PATH} "/home/jmwille/install/mpc-0.9/lib:/home/jmwille/install/mpfr-2.4.2/lib:/home/jmwille/install/gmp-4.3.2/lib:/home/jmwille/install/gcc-4.6.1/lib64:$ENV{LD_LIBRARY_PATH}")
SET(CTEST_TEST_TYPE EXPERIMENTAL)
#SET(CTEST_TEST_TIMEOUT 900)

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)
#Stokhos is explicitly disabled below to prevent the package from being
#implicitly enabled.  Sundance depends on Stokhos.
#SET(EXTRA_EXCLUDE_PACKAGES Phalanx Stokhos Sundance)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-DTPL_ENABLE_Pthread:BOOL=ON"
  "-DNOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK=ON"
  "-DMPI_BASE_DIR:PATH=/home/jmwille/install/gcc-4.6.1/openmpi-1.4.3"
  "-DBoost_INCLUDE_DIRS=/home/trilinos/tpl/gcc4.4.4/boost_1_46_1"
  "-DNetcdf_INCLUDE_DIRS=/home/trilinos/tpl/gcc4.4.4/netcdf-4.1.3/include"
  "-DNetcdf_LIBRARY_DIRS=/home/trilinos/tpl/gcc4.4.4/netcdf-4.1.3/lib"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

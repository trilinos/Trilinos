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


INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.trilinos-test.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME SERIAL_DEBUG_DEV2)
SET(COMPILER_VERSION "GCC-4.6.0")
SET(ENV{LD_LIBRARY_PATH} "/home/trilinos/install/gmp-4.3.2/lib:/home/trilinos/install/mpfr2.4.2/lib:/home/trilinos/install/mpc-0.8.1/lib:/home/trilinos/gcc4.6.0/lib64:$ENV{LD_LIBRARY_PATH}")
SET(CTEST_TEST_TYPE Nightly)
#SET(CTEST_TEST_TIMEOUT 900)

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)
#Disabling the following packages due to issues with them with the gc 4.6 compiler
#MOOCHO there were many test failures
#Piro there was one test failure
#Amesos there was one test failure
#Stratimikos there was one test failure
#Rythmos Disabling Stratimikos caused a lot of build errors
#stokhos there was one build error
SET(EXTRA_EXCLUDE_PACKAGES MOOCHO Piro Stratimikos Rythmos Zoltan2 Panzer)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DCMAKE_CXX_COMPILER:FILEPATH=/home/trilinos/gcc4.6.0/bin/g++"
  "-DCMAKE_C_COMPILER:FILEPATH=/home/trilinos/gcc4.6.0/bin/gcc"
  "-DCMAKE_Fortran_COMPILER:FILEPATH=/home/trilinos/gcc4.6.0/bin/gfortran"
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  "-DBoost_INCLUDE_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/boost_1_49_0"
  "-DTrilinos_ENABLE_MOOCHO:BOOL=OFF"
  "-DTrilinos_ENABLE_Piro:BOOL=OFF"
  "-DTrilinos_ENABLE_Stratimikos:BOOL=OFF"
  "-DTrilinos_ENABLE_Rythmos:BOOL=OFF"
  "-DHDF5_INCLUDE_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/hdf5-1.8.6/include"
  "-DHDF5_LIBRARY_DIRS:FILEPATH=/home/trilinos/tpl/gcc4.1.2/hdf5-1.8.6/lib"
  "-DNetcdf_LIBRARY_DIRS=/home/trilinos/tpl/gcc4.1.2/netcdf_4.2/lib"
  "-DNetcdf_INCLUDE_DIRS=/home/trilinos/tpl/gcc4.1.2/netcdf_4.2/include"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

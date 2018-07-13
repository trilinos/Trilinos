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


INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.s956330.clang3.4.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME SERIAL_OPT_DEV_CLANG)
SET(COMPILER_VERSION "Clang-3.4")
SET(ENV{LD_LIBRARY_PATH} "/projects/sems/compilers/clang/3.4/lib:$ENV{LD_LIBRARY_PATH}")
SET(CTEST_TEST_TYPE Experimental)
#SET(CTEST_TEST_TIMEOUT 900)
SET( CTEST_BUILD_FLAGS "-j11 -i" )
SET( CTEST_PARALLEL_LEVEL "11" )

SET(Trilinos_ENABLE_SECONDARY_TESTED_CODE ON)
#Stokhos is explicitly disabled below to prevent the package from being
#implicitly enabled.  Sundance depends on Stokhos.
SET(EXTRA_EXCLUDE_PACKAGES PyTrilinos Optika )

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-DTPL_ENABLE_Pthread:BOOL=ON"
  "-DTPL_ENABLE_Boost:BOOL=ON"
  "-DBoost_INCLUDE_DIRS:FILEPATH=/projects/sems/tpls/gcc_4.8.2/arch_x86/boost_1.49.0"
  "-DTPL_ENABLE_SuperLU=ON"
  "-DSuperLU_INCLUDE_DIRS:PATH=/projects/sems/tpls/clang_3.4/arch_x86_64/SuperLU_4.3/SRC"
  "-DSuperLU_LIBRARY_DIRS:PATH=/projects/sems/tpls/clang_3.4/arch_x86_64/SuperLU_4.3/lib"
  "-DSuperLU_LIBRARY_NAMES:STRING=superlu_4.3"
  "-DNOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK=ON"
  "-DCMAKE_CXX_COMPILER:FILEPATH=/projects/sems/compilers/clang/3.4/bin/clang++"
  "-DCMAKE_C_COMPILER:FILEPATH=/projects/sems/compilers/clang/3.4/bin/clang"
  "-DCMAKE_Fortran_COMPILER:FILEPATH=/projects/sems/compilers/gcc/4.8.2/bin/gfortran"
  "-DTPL_ENABLE_HDF5:BOOL=OFF"
  "-DHDF5_INCLUDE_DIRS:FILEPATH=/projects/sems/tpls/gcc_4.8.2/arch_x86/hdf5_1.8.12/include"
  "-DHDF5_LIBRARY_DIRS:FILEPATH=/projects/sems/tpls/gcc_4.8.2/arch_x86/hdf5_1.8.12/lib"
  "-DTPL_ENABLE_Netcdf:BOOL=OFF"
  "-DNetcdf_LIBRARY_DIRS=/projects/sems/tpls/gcc_4.8.2/arch_x86/pnetcdf_4.3.0/lib"
  "-DNetcdf_INCLUDE_DIRS=/projects/sems/tpls/gcc_4.8.2/arch_x86/pnetcdf_4.3.0/include"
  "-DSWIG_EXECUTABLE:FILEPATH=/projects/sems/tpls/gcc_4.8.2/arch_x86/swig_2.0.11/bin/swig"
  "-DTPL_ENABLE_METIS:BOOL=OFF"
  "-DMETIS_LIBRARY_DIRS:PATH=/projects/sems/tpls/gcc_4.8.2/arch_x86/metis_5.1.0/lib"
  "-DMETIS_INCLUDE_DIRS:PATH=/projects/sems/tpls/gcc_4.8.2/arch_x86/metis_5.1.0/include"
  "-DZoltan_ENABLE_ParMETIS=OFF"
  "-DZoltan_ENABLE_METIS=OFF"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()

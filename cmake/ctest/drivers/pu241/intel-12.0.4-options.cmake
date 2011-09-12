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

INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/mkl-12.0.4-options.cmake)
SET(BLAS_INCLUDE_DIRS   ${MKLROOT}/include/intel64/lp64                     CACHE PATH  "Path to MKL BLAS Fortran modules compatible with Intel fortran")
SET(LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS}                                CACHE PATH  "Path to MKL LAPACK Fortran modules compatible with Intel fortran")

INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/boost-1.46.1-options.cmake)
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/casl-vri-tpls.cmake)

SET(PVMLibraries_LIBRARY_DIRS /opt/intel-11.1.064/tpls/pvm3/lib/LINUX64                    CACHE FILEPATH "")
SET(PVMHeaders_INCLUDE_DIRS /opt/intel-11.1.064/tpls/pvm3/include                        CACHE FILEPATH "")
SET(HDF5_LIBRARY_NAMES "hdf5_hl;hdf5;hdf5_cpp"                                    CACHE STRING   "")
SET(HDF5_LIBRARY_DIRS /opt/intel-11.1.064/tpls/hdf5-1.8.5-patch1/lib              CACHE FILEPATH "")
SET(HDF5_INCLUDE_DIRS /opt/intel-11.1.064/tpls/hdf5-1.8.5-patch1/include          CACHE FILEPATH "")
SET(EpetraExt_ENABLE_HDF5 OFF CACHE BOOL "")

SET(INTEL_BIN /opt/intel/Compiler/composerxe-2011.4.191/bin/intel64)

SET(CMAKE_SKIP_RPATH ON BOOL "")
#SET(BUILD_SHARED_LIBS ON CACHE BOOL "")
SET(TPL_ENABLE_BinUtils ON CACHE BOOL "")
SET(CMAKE_C_COMPILER ${INTEL_BIN}/icc CACHE FILEPATH "")
SET(CMAKE_CXX_COMPILER ${INTEL_BIN}/icpc CACHE FILEPATH "")
SET(CMAKE_Fortran_COMPILER ${INTEL_BIN}/ifort CACHE FILEPATH "")

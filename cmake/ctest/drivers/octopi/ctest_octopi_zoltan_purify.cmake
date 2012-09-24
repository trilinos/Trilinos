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


INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")


INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../../../Version.cmake")

#
# Set the options specific to this build case
#

SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME "ZOLTAN_PURIFY")
SET(Trilinos_TRACK ${Trilinos_TESTING_TRACK})
SET(Trilinos_BRANCH ${Trilinos_REPOSITORY_BRANCH})
SET(CTEST_TEST_TYPE Nightly)
SET(CTEST_TEST_TIMEOUT "2400")

#
# Set the rest of the system-specific options and run the dashboard build/test
# Platform/compiler specific options for octopi
#

# Base of Trilinos/cmake/ctest then BUILD_DIR_NAME

SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )

SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )

SET( CTEST_BUILD_FLAGS "-j8 -i" )

SET_DEFAULT( CTEST_PARALLEL_LEVEL "1" )
SET_DEFAULT(COMPILER_VERSION "GCC-3.4.6")


SET_DEFAULT( Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)

# Only turn on PyTrilinos for shared libraries
SET_DEFAULT(Trilinos_EXCLUDE_PACKAGES TriKota Optika)

# Purify info
SET(KDD_PURIFY "/usr/local/rational/releases/PurifyPlus.7.0/i386_linux2/bin/purify")
SET(KDD_PURIFY_ARGS "-best-effort -follow-child-processes=yes -cache-dir=/tmp/purify -chain-length=20  -windows=no ")

# Compilers that work with purify
SET(KDD_GCC  "/usr/bin/gcc346")
SET(KDD_GCXX "/usr/bin/g++346")

# Output of "mpicc --showme:compile" and "mpiCC --showme:compile"
SET(KDD_CFLAGS   "-std=c99 -m64 -g -I/opt/lam714-gcc346-pure/include -pthread")
SET(KDD_CXXFLAGS "-m64 -g -I/opt/lam714-gcc346-pure/include -pthread")

# Output of "mpiCC --showme:link"
SET(KDD_LINKFLAGS "-m64 -L/opt/lam714-gcc346-pure/lib -llammpio -llammpi++ -llamf77mpi -lmpi -llam -laio -laio -lutil -ldl")

# MPI info; needed for mpirun; also need this in path.
set(KDD_MPI "/opt/lam714-gcc346-pure/bin")
set(TPL_ENABLE_MPI ON)

set(TDD_HTTP_PROXY "http://sonproxy.sandia.gov:80/")

SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
    "-DTPL_ENABLE_MPI:BOOL=ON"
    "-DMPI_USE_COMPILER_WRAPPERS:BOOL=OFF"
    "-DMPI_BIN_DIR:STRING=${KDD_MPI}"
    "-DCMAKE_C_COMPILER:STRING=${KDD_PURIFY}" 
    "-DCMAKE_C_FLAGS:STRING=${KDD_PURIFY_ARGS} ${KDD_GCC} ${KDD_CFLAGS}"
    "-DCMAKE_CXX_COMPILER:STRING=${KDD_PURIFY}" 
    "-DCMAKE_CXX_FLAGS:STRING=${KDD_PURIFY_ARGS} ${KDD_GCXX} ${KDD_CXXFLAGS}"
    "-DTrilinos_EXTRA_LINK_FLAGS:STRING=${KDD_LINKFLAGS}"
    "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON "
    "-DMPI_EXEC_MAX_NUMPROCS:STRING=11 "
    "-DTrilinos_ENABLE_Fortran:BOOL=OFF "
    "-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF "
    "-DTrilinos_ENABLE_EXAMPLES:BOOL=ON "
    "-DTrilinos_VERBOSE_CONFIGURE:BOOL=ON "
    "-DTrilinos_ENABLE_Zoltan:BOOL=ON "
    "-DZoltan_ENABLE_EXAMPLES:BOOL=OFF "
    "-DZoltan_ENABLE_TESTS:BOOL=ON "
    "-DZoltan_ENABLE_ParMETIS:BOOL=ON "
    "-DParMETIS_LIBRARY_DIRS:FILEPATH=/Net/local/proj/zoltan/arch/linux64/lib/lam/ParMETIS3" 
    "-DParMETIS_INCLUDE_DIRS:FILEPATH=/Net/local/proj/zoltan/arch/all/src/ParMETIS3" 
    "-DZoltan_ENABLE_Scotch:BOOL=ON"
    "-DScotch_LIBRARY_DIRS:FILEPATH=/Net/local/proj/zoltan/arch/linux64/lib/lam/Scotch5" 
    "-DScotch_INCLUDE_DIRS:FILEPATH=/Net/local/proj/zoltan/arch/all/src/Scotch5" 
    "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
    "-DTrilinos_ENABLE_DEPENCENCY_UNIT_TESTS:BOOL=OFF"
  )

TRILINOS_CTEST_DRIVER()

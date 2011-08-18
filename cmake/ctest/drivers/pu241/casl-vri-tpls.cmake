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

SET(TPL_ENABLE_JDK     ON      CACHE BOOL "")
SET(TPL_ENABLE_OpenSSL ON      CACHE BOOL "")
SET(TPL_ENABLE_Zlib    ON      CACHE BOOL "")
SET(TPL_ENABLE_TCL     ON      CACHE BOOL "")
SET(TPL_ENABLE_Netcdf ON CACHE BOOL "")

# We don't have the MATLAB TPL for SEACAS
SET(TPL_ENABLE_MATLAB OFF CACHE BOOL "" FORCE)

SET(ANT_PATH         /usr/bin                                                     CACHE FILEPATH "")
SET(TCLSH_PATH       /usr/bin                                                     CACHE FILEPATH "")
SET(TCL_LIBRARY_NAMES "tcl8.5"                                                    CACHE STRING   "")
SET(JDK_INCLUDE_DIRS "/usr/lib/jvm/java/include;/usr/lib/jvm/java/include/linux"  CACHE FILEPATH "")
SET(JDK_LIBRARY_DIRS /usr/lib/jvm/java/jre/lib/amd64/server                       CACHE FILEPATH "")
SET(Netcdf_INCLUDE_DIRS /opt/gcc-4.5.1/tpls/netcdf-4.1.1/include CACHE FILEPATH "")
SET(Netcdf_LIBRARY_DIRS /opt/gcc-4.5.1/tpls/netcdf-4.1.1/lib CACHE FILEPATH "")
SET(SILO_INCLUDE_DIRS   /opt/gcc-4.5.1/tpls/silo/include CACHE FILEPATH "")
SET(SILO_LIBRARY_DIRS   /opt/gcc-4.5.1/tpls/silo/lib     CACHE FILEPATH "")
SET(TBB_INCLUDE_DIRS    /opt/intel/Compiler/composerxe-2011.4.191/tbb/include                                      CACHE FILEPATH "")
SET(TBB_LIBRARY_DIRS    /opt/intel/Compiler/composerxe-2011.4.191/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21  CACHE FILEPATH "")
SET(LIBXML2_INCLUDE_DIRS  /usr/include/libxml2       CACHE FILEPATH "")
SET(LIBXML2_LIBRARY_DIRS  /usr/lib64                  CACHE FILEPATH "")

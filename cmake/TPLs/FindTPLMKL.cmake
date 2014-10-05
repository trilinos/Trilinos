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

INCLUDE(TribitsTplDeclareLibraries)

#-----------------------------------------------------------------------------
#  Intel's Math Kernel Library (MKL)
#
#  Acquisition information:
#    Date checked:  06 Aug 2012
#    Checked by:    Mark Hoemmen <mhoemme AT sandia.gov>
#    Version:       10.3
#

# Intel's Math Kernel Library (MKL) provides an implementation of the
# BLAS and LAPACK, which can be used to satisfy the BLAS and LAPACK
# TPLs.  However, the "MKL TPL" here refers to other functionality
# provided by the MKL, including things like sparse matrix kernels and
# pseudorandom number generators.  That's why we require a header
# file, to access the function declarations.

TRIBITS_TPL_DECLARE_LIBRARIES( MKL
  REQUIRED_HEADERS mkl.h
  REQUIRED_LIBS_NAMES mkl_rt
  )

# In the past, MKL users had to link with a long list of libraries.
# The choice of libraries enables specific functionality.  Intel
# provides a web page (the "Intel MKL Link Line Advisor") to help
# users pick which libraries and compiler flags to use:
#
# http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/
#
# As of version 10.3, MKL has an option to use a single dynamic
# library via "-lmkl_rt".  See the following article:
#
# http://software.intel.com/en-us/articles/a-new-linking-model-single-dynamic-library-mkl_rt-since-intel-mkl-103/
#
# This is why we made mkl_rt the required library name.  Users must
# override this if they are using an older MKL version, or if they
# prefer the long list of libraries to the single dynamic library
# option.  (I highly recommend the single dynamic library option, if
# your system allows dynamic libraries.)
#
# Users will probably need to specify the CMake options
# MKL_LIBRARY_DIRS (the path to the libraries) and MKL_INCLUDE_DIRS
# (the path to the header files).  On Linux, you may also have to link
# with Pthreads and libm (the C math library).  This may happen by
# default, but if you have trouble linking, try setting
# MKL_LIBRARY_NAMES to "mkl_rt;pthread;m" (for the single dynamic
# library option).  If you still have trouble, or if you are unable to
# use the single dynamic library option, look at examples of linking
# with MKL as the BLAS and LAPACK implementation in the sampleScripts/
# subdirectory of the Trilinos source tree.  Copy BLAS_LIBRARY_NAMES
# to MKL_LIBRARY_NAMES and try again.
#
# Users may also need to use special compiler flags, in particular if
# they wish to enable ILP64 support (for 64-bit integer indices).  See
# the above "Intel MKL Link Line Advisor" link for details.
#
# Here is a combination of CMake options that worked for me, when
# building and running on Linux, using the single dynamic library
# option with MKL 10.3:
#
# -D BLAS_LIBRARY_DIRS:FILEPATH="${MKLROOT}/lib/intel64" 
# -D BLAS_LIBRARY_NAMES:STRING="mkl_rt" 
# -D LAPACK_LIBRARY_DIRS:FILEPATH="${MKLROOT}/lib/intel64" 
# -D LAPACK_LIBRARY_NAMES:STRING="mkl_rt" 
# -D MKL_LIBRARY_DIRS:FILEPATH="${MKLROOT}/lib/intel64" 
# -D MKL_LIBRARY_NAMES:STRING="mkl_rt" 
# -D MKL_INCLUDE_DIRS:FILEPATH="${MKLROOT}/include" 
# 
# where the MKLROOT environment variable points to my MKL install
# directory.






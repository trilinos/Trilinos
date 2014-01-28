/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_MACROS_HPP
#define KOKKOS_MACROS_HPP

#include <KokkosCore_config.h>

namespace Kokkos {
class HostSpace ;
class CudaSpace ;
} // namespace Kokkos

//----------------------------------------------------------------------------
//-------- Identify Compiler Version -----------------------------------------
//----------------------------------------------------------------------------

#if defined __ECC || defined __ICC || defined __INTEL_COMPILER
  #define KOKKOS_COMPILER_NAME "Intel C++"
  #if defined __ICC
    #define KOKKOS_COMPILER_VERSION __ICC
  #else
    #if defined __INTEL_COMPILER
      #define KOKKOS_COMPILER_VERSION __INTEL_COMPILER
    #else
      #define KOKKOS_COMPILER_VERSION __ECC
    #endif
  #endif
  #define KOKKOS_COMPILER_INTEL 1
#endif

#if defined __IBMC__ || defined __IBMCPP__
  #define KOKKOS_COMPILER_NAME "IBM C++"
  #if defined __IBMC__
    #define KOKKOS_COMPILER_VERSION __IBMC__
  #else
    #define KOKKOS_COMPILER_VERSION __IBMCPP__
  #endif
  #define KOKKOS_COMPILER_IBM 1
#endif

#if defined __APPLE_CC__
   /* Apple uses GNU as compiler */
  #define KOKKOS_COMPILER_APPLECC 1
#endif

#if defined __clang__
  #define KOKKOS_COMPILER_NAME "Clang"
  #define KOKKOS_COMPILER_VERSION __clang_major__*100+__clang_minor__*10+__clang_patchlevel__
  #define KOKKOS_COMPILER_CLANG 1
#endif

#if defined __GNUC__ && !defined KOKKOS_COMPILER_NAME && !defined __clang__
  #define KOKKOS_COMPILER_NAME "Gnu GCC"
  #define KOKKOS_COMPILER_VERSION __GNUC__*100+__GNUC_MINOR__*10+__GNUC_PATCHLEVEL__
  #define KOKKOS_COMPILER_GCC 1
#endif

#if defined __PGIC__ && !defined KOKKOS_COMPILER_NAME
  #define KOKKOS_COMPILER_NAME "PGI C++"
  #define KOKKOS_COMPILER_VERSION __PGIC__*100+__PGIC_MINOR__*10+__PGIC_PATCHLEVEL__
  #define KOKKOS_COMPILER_PGI 1
#endif

#if defined __NVCC__
  #define KOKKOS_DEVICE_COMPILER_NAME "NVIDIA NVCC"
  #define KOKKOS_DEVICE_COMPILER_VERSION __NVCC__
#endif

#if !defined KOKKOS_COMPILER_NAME
  #define KOKKOS_COMPILER_NAME "Unknown compiler"
#endif

#if !defined KOKKOS_COMPILER_VERSION
  #define KOKKOS_COMPILER_VERSION 0
#endif

#if !defined KOKKOS_DEVICE_COMPILER_NAME
  #define KOKKOS_DEVICE_COMPILER_NAME KOKKOS_COMPILER_NAME
#endif

#if !defined KOKKOS_DEVICE_COMPILER_VERSION
  #define KOKKOS_DEVICE_COMPILER_VERSION KOKKOS_COMPILER_VERSION
#endif


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( __CUDACC__ ) 

// Compiling with CUDA compiler.

#if ! defined( KOKKOS_HAVE_CUDA )
#error "Compiling Kokkos with Cuda compiler but KOKKOS_HAVE_CUDA is undefined"
#endif

#include <cuda.h>

/*  Compiling with a CUDA compiler for device code.
 *
 *  Include <cuda.h> to pick up the CUDA_VERSION macro defined as:
 *    CUDA_VERSION = ( MAJOR_VERSION * 1000 ) + ( MINOR_VERSION * 10 )
 *
 *  When generating device code the __CUDA_ARCH__ macro is defined as:
 *    __CUDA_ARCH__ = ( MAJOR_CAPABILITY * 100 ) + ( MINOR_CAPABILITY * 10 )
 */
#if ! defined( CUDA_VERSION )
#error "#include <cuda.h> did not define CUDA_VERSION"
#endif

#if ( CUDA_VERSION < 4010 )
#error "Cuda version 4.1 or greater required"
#endif

#endif /* #if defined( __CUDACC__ ) */

//----------------------------------------------------------------------------

#if defined( __CUDACC__ ) && defined( __CUDA_ARCH__ )

/*  Compiling with CUDA compiler for device code. */

#if ( __CUDA_ARCH__ < 200 )
#error "Cuda device capability >= 2.0 is required"
#endif

#define KOKKOS_FORCEINLINE_FUNCTION  __device__  __host__  __forceinline__
#define KOKKOS_INLINE_FUNCTION  __device__  __host__  inline
#define KOKKOS_FUNCTION         __device__  __host__

#endif /* #if defined( __CUDACC__ ) && #if defined( __CUDA_ARCH__ ) */

//----------------------------------------------------------------------------

#if defined( __CUDACC__ ) && ! defined( __CUDA_ARCH__ )

/*  Compiling with CUDA compiler for host code. */

#define KOKKOS_FORCEINLINE_FUNCTION  __forceinline__

#endif /* #if defined( __CUDACC__ ) && ! defined( __CUDA_ARCH__ ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( __INTEL_COMPILER )

#if (__INTEL_COMPILER < 1200)
#define KOKKOS_DISABLE_ASM true;
#endif

/*  Compiling with Intel compiler */
/*  TBD: Version testing */

#ifndef KOKKOS_FORCEINLINE_FUNCTION
#define KOKKOS_FORCEINLINE_FUNCTION  __forceinline
#endif

#if defined( __MIC__ )

/*  Compiling with Intel compiler for execution on an Intel MIC device.
 *  These devices are used in no-offload mode so the HostSpace is the MIC space.
 */

#else

#ifndef KOKKOS_USE_PRAGMA_SIMD
#define KOKKOS_USE_PRAGMA_SIMD
#endif

/*
  #pragma simd vectorlength(N)
  #pragma ivdep
*/

#endif /* #if defined( __MIC__ ) */

#endif /* #if defined( __INTEL_COMPILER ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( __GNUC__ ) /* GNU C   */ || \
    defined( __GNUG__ ) /* GNU C++ */

/* Compiling with GNU compiler */

#ifndef KOKKOS_FORCEINLINE_FUNCTION
#define KOKKOS_FORCEINLINE_FUNCTION  inline __attribute__((always_inline))
#endif

/*  Compiling with GNU compatible compiler.  */

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( _OPENMP )

#if ! defined( KOKKOS_HAVE_OPENMP )
#error "Compiling Kokkos for OpenMP but KOKKOS_HAVE_OPENMP is undefined"
#endif

/*  Compiling with OpenMP.
 *  The value of _OPENMP is an integer value YYYYMM
 *  where YYYY and MM are the year and month designation
 *  of the supported OpenMP API version.
 */

#endif /* END: #if defined( _OPENMP ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef KOKKOS_FUNCTION
#define KOKKOS_FUNCTION /* */
#endif

#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif

#ifndef KOKKOS_FORCEINLINE_FUNCTION
#define KOKKOS_FORCEINLINE_FUNCTION  inline
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( __CUDACC__ ) && defined( __CUDA_ARCH__ )

namespace Kokkos { typedef CudaSpace ExecutionSpace ; }

#else

namespace Kokkos { typedef HostSpace ExecutionSpace ; }

#endif

#define KOKKOS_RESTRICT_EXECUTION_TO_DATA( DATA_SPACE , DATA_PTR ) \
  Kokkos::VerifyExecutionSpaceCanAccessDataSpace< \
    Kokkos::ExecutionSpace , DATA_SPACE >::verify( DATA_PTR )

#define KOKKOS_RESTRICT_EXECUTION_TO( DATA_SPACE ) \
  Kokkos::VerifyExecutionSpaceCanAccessDataSpace< \
    Kokkos::ExecutionSpace , DATA_SPACE >::verify()

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_MACROS_HPP */


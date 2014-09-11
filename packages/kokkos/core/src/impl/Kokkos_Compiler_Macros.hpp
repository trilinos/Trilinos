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

#ifndef KOKKOS_COMPILER_MACROS_HPP
#define KOKKOS_COMPILER_MACROS_HPP

#if defined( KOKKOS_HAVE_CUDA )
//----------------------------------------------------------------------------

#if defined( __CUDACC__ )

/*  Compiling with a CUDA compiler.
 *
 *  Include <cuda.h> to pick up the CUDA_VERSION macro defined as:
 *    CUDA_VERSION = ( MAJOR_VERSION * 1000 ) + ( MINOR_VERSION * 10 )
 *
 *  When generating device code the __CUDA_ARCH__ macro is defined as:
 *    __CUDA_ARCH__ = ( MAJOR_CAPABILITY * 100 ) + ( MINOR_CAPABILITY * 10 )
 */

#include <cuda_runtime.h>
#include <cuda.h>

#if ! defined( CUDA_VERSION )
#error "#include <cuda.h> did not define CUDA_VERSION"
#endif

#if ( CUDA_VERSION < 4010 )
#error "Cuda version 4.1 or greater required"
#endif

#if defined( __CUDA_ARCH__ ) && ( __CUDA_ARCH__ < 200 )
/*  Compiling with CUDA compiler for device code. */
#error "Cuda device capability >= 2.0 is required"
#endif

#endif /* #if defined( __CUDACC__ ) */

//----------------------------------------------------------------------------

#endif

/*--------------------------------------------------------------------------*/
/* Language info: C++, CUDA, OPENMP */

#if defined( __CUDA_ARCH__ )
  // Compiling Cuda code to 'ptx'

  #define KOKKOS_FORCEINLINE_FUNCTION  __device__  __host__  __forceinline__
  #define KOKKOS_INLINE_FUNCTION       __device__  __host__  inline
  #define KOKKOS_FUNCTION              __device__  __host__

#endif /* #if defined( __CUDA_ARCH__ ) */

#if defined( _OPENMP )

  /*  Compiling with OpenMP.
   *  The value of _OPENMP is an integer value YYYYMM
   *  where YYYY and MM are the year and month designation
   *  of the supported OpenMP API version.
   */

#endif /* #if defined( _OPENMP ) */

/*--------------------------------------------------------------------------*/
/* Mapping compiler built-ins to KOKKOS_COMPILER_*** macros */

#if defined( __NVCC__ )
  // NVIDIA compiler is being used.
  // Code is parsed and separated into host and device code.
  // Host code is compiled again with another compiler.
  // Device code is compile to 'ptx'.
  #define KOKKOS_COMPILER_NVCC __NVCC__

  #if defined( KOKKOS_HAVE_CXX11 )
    // CUDA supports (inofficially) C++11 in device code starting with 
    // version 6.5. This includes auto type and device code internal
    // lambdas.
    #if ( CUDA_VERSION < 6050 )
      #error "NVCC does not support C++11"
    #endif
  #endif
#else
  #if defined( KOKKOS_HAVE_CXX11 )
    // CUDA (including version 6.5) does not support giving lambdas as
    // arguments to global functions. Thus its not currently possible
    // to dispatch lambdas from the host.
    #define KOKKOS_HAVE_CXX11_DISPATCH_LAMBDA 1
  #endif
#endif /* #if defined( __NVCC__ ) */


#if ! defined( __CUDA_ARCH__ ) /* Not compiling Cuda code to 'ptx'. */

#if defined( __INTEL_COMPILER )
  #define KOKKOS_COMPILER_INTEL __INTEL_COMPILER
#elif defined( __ICC )
  // Old define
  #define KOKKOS_COMPILER_INTEL __ICC
#elif defined( __ECC ) 
  // Very old define
  #define KOKKOS_COMPILER_INTEL __ECC
#endif

#if defined( _CRAYC )
  #define KOKKOS_COMPILER_CRAYC _CRAYC
#endif

#if defined( __IBMCPP__ )
  // IBM C++
  #define KOKKOS_COMPILER_IBM __IBMCPP__
#elif defined( __IBMC__ )
  #define KOKKOS_COMPILER_IBM __IBMC__
#endif

#if defined( __APPLE_CC__ )
  #define KOKKOS_COMPILER_APPLECC __APPLE_CC__
#endif

#if defined( __clang__ )
  #define KOKKOS_COMPILER_CLANG __clang_major__*100+__clang_minor__*10+__clang_patchlevel__
#endif

#if ! defined( __clang__ ) && ! defined( KOKKOS_COMPILER_INTEL ) &&defined( __GNUC__ )
  #define KOKKOS_COMPILER_GNU __GNUC__*100+__GNUC_MINOR__*10+__GNUC_PATCHLEVEL__
#endif

#if defined( __PGIC__ ) && ! defined( __GNUC__ )
  #define KOKKOS_COMPILER_PGI __PGIC__*100+__PGIC_MINOR__*10+__PGIC_PATCHLEVEL__
#endif

#endif /* #if ! defined( __CUDA_ARCH__ ) */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/* Intel compiler macros */

#if defined( KOKKOS_COMPILER_INTEL )

  #define KOKKOS_HAVE_PRAGMA_UNROLL 1
  #define KOKKOS_HAVE_PRAGMA_IVDEP 1
  #define KOKKOS_HAVE_PRAGMA_LOOPCOUNT 1
  #define KOKKOS_HAVE_PRAGMA_VECTOR 1
  #define KOKKOS_HAVE_PRAGMA_SIMD 1

  #if ( 1200 <= KOKKOS_COMPILER_INTEL ) && ! defined( KOKKOS_ENABLE_ASM )
    #define KOKKOS_ENABLE_ASM 1
  #endif

  #define KOKKOS_FORCEINLINE_FUNCTION  __forceinline

  #if defined( __MIC__ )
    // Compiling for Xeon Phi
  #endif

#endif

/*--------------------------------------------------------------------------*/
/* Cray compiler macros */

#if defined( KOKKOS_COMPILER_CRAYC )


#endif

/*--------------------------------------------------------------------------*/
/* IBM Compiler macros */

#if defined( KOKKOS_COMPILER_IBM )

  #define KOKKOS_HAVE_PRAGMA_UNROLL 1
  //#define KOKKOS_HAVE_PRAGMA_IVDEP 1
  //#define KOKKOS_HAVE_PRAGMA_LOOPCOUNT 1
  //#define KOKKOS_HAVE_PRAGMA_VECTOR 1
  //#define KOKKOS_HAVE_PRAGMA_SIMD 1

#endif

/*--------------------------------------------------------------------------*/

#if defined( KOKKOS_COMPILER_CLANG )

  //#define KOKKOS_HAVE_PRAGMA_UNROLL 1
  //#define KOKKOS_HAVE_PRAGMA_IVDEP 1
  //#define KOKKOS_HAVE_PRAGMA_LOOPCOUNT 1
  //#define KOKKOS_HAVE_PRAGMA_VECTOR 1
  //#define KOKKOS_HAVE_PRAGMA_SIMD 1

  #define KOKKOS_FORCEINLINE_FUNCTION  inline __attribute__((always_inline))

#endif

/*--------------------------------------------------------------------------*/

#if defined( KOKKOS_COMPILER_GNU ) 

  //#define KOKKOS_HAVE_PRAGMA_UNROLL 1
  //#define KOKKOS_HAVE_PRAGMA_IVDEP 1
  //#define KOKKOS_HAVE_PRAGMA_LOOPCOUNT 1
  //#define KOKKOS_HAVE_PRAGMA_VECTOR 1
  //#define KOKKOS_HAVE_PRAGMA_SIMD 1

  #define KOKKOS_FORCEINLINE_FUNCTION inline __attribute__((always_inline))

  #if ! defined( KOKKOS_ENABLE_ASM ) && \
      ! ( defined( __powerpc) || \
          defined(__powerpc__) || \
          defined(__powerpc64__) || \
          defined(__POWERPC__) || \
          defined(__ppc__) || \
          defined(__ppc64__) )
    #define KOKKOS_ENABLE_ASM 1
  #endif

  #define KOKKOS_NONTEMPORAL_PREFETCH_LOAD(addr) __builtin_prefetch(addr,0,0)
  #define KOKKOS_NONTEMPORAL_PREFETCH_STORE(addr) __builtin_prefetch(addr,1,0)

#endif

/*--------------------------------------------------------------------------*/

#if defined( KOKKOS_COMPILER_PGI )

  #define KOKKOS_HAVE_PRAGMA_UNROLL 1
  #define KOKKOS_HAVE_PRAGMA_IVDEP 1
  //#define KOKKOS_HAVE_PRAGMA_LOOPCOUNT 1
  #define KOKKOS_HAVE_PRAGMA_VECTOR 1
  //#define KOKKOS_HAVE_PRAGMA_SIMD 1

#endif

/*--------------------------------------------------------------------------*/

#if defined( KOKKOS_COMPILER_NVCC )

  #if defined(__CUDA_ARCH__ )
    #define KOKKOS_HAVE_PRAGMA_UNROLL 1
  #endif

#endif

/*--------------------------------------------------------------------------*/
/* Select compiler dependent interface for atomics */

#if ! defined( KOKKOS_ATOMICS_USE_CUDA ) || \
    ! defined( KOKKOS_ATOMICS_USE_GNU ) || \
    ! defined( KOKKOS_ATOMICS_USE_INTEL ) || \
    ! defined( KOKKOS_ATOMICS_USE_OPENMP31 )

/* Atomic selection is not pre-defined, choose from language and compiler. */

#if defined( __CUDA_ARCH__ )

  #define KOKKOS_ATOMICS_USE_CUDA

#elif defined( KOKKOS_COMPILER_GNU ) || defined( KOKKOS_COMPILER_CLANG )

  #define KOKKOS_ATOMICS_USE_GNU

#elif defined( KOKKOS_COMPILER_INTEL ) || defined( KOKKOS_COMPILER_CRAYC )

  #define KOKKOS_ATOMICS_USE_INTEL

#elif defined( _OPENMP ) && ( 201107 <= _OPENMP )

  #define KOKKOS_ATOMICS_USE_OMP31

#else

  #error "Compiler does not support atomic operations"

#endif

#endif

/*--------------------------------------------------------------------------*/

#endif /* #ifndef KOKKOS_COMPILER_MACROS_HPP */


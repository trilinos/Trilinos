// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef CONFIG_STK_SIMD_H
#define CONFIG_STK_SIMD_H

// IWYU pragma: private, include <stk_simd/Simd.hpp>

#if !defined(__CUDA_ARCH__)
//#if 0 // to turn off simd

#if ( defined(_M_AMD64) || defined(_M_X64) || defined(__amd64) ) && ! defined(__x86_64__)
#define __x86_64__
#endif

#if defined(__SSE2__) || defined(__SSE3__) || defined(__SSSE3__) || (_M_IX86_FP >= 2)
#define __SSE_23__
#endif

#if defined(__SSE4_1__) || defined(__SSE4_2__)
#define __SSE4_X__
#endif

#if defined(__AVX__) || defined(__AVX2__)
#define __AVX_X__
#endif

#if defined(__AVX512F__)
#define __AVX512_X__
#endif

#else // using cuda

#ifndef USE_STK_SIMD_NONE
#define USE_STK_SIMD_NONE
#endif

#endif

/* cpp defines are processed with the following priority:
   1. USE_STK_SIMD_NONE
   2. USE_STK_SIMD_SSE
   3. USE_STK_SIMD_AVX
   4. USE_STK_SIMD_AUTO
*/

// default to using AUTO

#if defined (USE_STK_SIMD_NONE) || (USE_STK_SIMD_SSE) || (USE_STK_SIMD_AVX) || (USE_STK_SIMD_AVX512)
#else
#undef  USE_STK_SIMD_AUTO
#define USE_STK_SIMD_AUTO // default to SSE or less capability
#endif

#ifdef USE_STK_SIMD_NONE

#  define STK_SIMD_NONE

#else

#  if defined (USE_STK_SIMD_SSE)
#    if defined(__SSE_23__) || defined(__SSE4_X__) || defined(__x86_64__)
#      define STK_SIMD_SSE
#    else
#      error "STK: Build asked for SIMD/SSE but processor does not support it"
#    endif
#  elif defined (USE_STK_SIMD_AVX)
#    if defined(__AVX_X__)
#      define STK_SIMD_AVX
#    else
#      error "STK: Build asked for SIMD/AVX but processor does not support it"
#    endif
#  elif defined (USE_STK_SIMD_AVX512)
#    define STK_SIMD_AVX512
#  elif defined (USE_STK_SIMD_AUTO)
#    if defined(__AVX512_X__)
#      define STK_SIMD_AVX512
#    elif defined(__AVX_X__)
#      define STK_SIMD_AVX
#    elif defined(__SSE_23__) || defined(__SSE4_X__) || defined(__x86_64__)
#      define STK_SIMD_SSE
#    elif defined(__ARM_NEON) || defined(__VSX__)
#    else
#      define STK_SIMD_NONE
#    endif
#  endif

#  if defined(STK_SIMD_SSE) && defined(STK_SIMD_AVX)
#    error "STK: Cannot have both STK_SIMD_SSE and STK_SIMD_AVX defined"
#  endif

#  if defined(STK_SIMD_SSE) && defined(STK_SIMD_AVX512)
#    error "STK: Cannot have both STK_SIMD_SSE and STK_SIMD_AVX512 defined"
#  endif

#  if defined(STK_SIMD_AVX) && defined(STK_SIMD_AVX512)
#    error "STK: Cannot have both STK_SIMD_AVX and STK_SIMD_AVX512 defined"
#  endif

#endif

#if defined (STK_KOKKOS_SIMD)
#undef STK_SIMD_SSE
#undef STK_SIMD_AVX
#undef STK_SIMD_AVX512
#endif

#endif // #ifndef CONFIG_STK_SIMD_H


#ifndef CONFIG_STK_SIMD_H
#define CONFIG_STK_SIMD_H

#if !defined(__CUDA_ARCH__)

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

//## #if defined(__AVX512F__)
//## #define __AVX512_X__
//## #endif

#else // using cuda

#define USE_STK_SIMD_NONE

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

#  define STK_NO_SIMD
#  ifdef STK_SIMD
#    undef STK_SIMD
#  endif

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
#    endif
#  endif

#  if defined(STK_SIMD_SSE) || defined(STK_SIMD_AVX) || defined(STK_SIMD_AVX512)
#    define STK_SIMD
#    ifdef STK_NO_SIMD
#      undef STK_NO_SIMD
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

#endif // #ifndef CONFIG_STK_SIMD_H






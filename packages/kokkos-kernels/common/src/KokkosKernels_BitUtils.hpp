//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef _KOKKOSKERNELS_BITUTILS_HPP
#define _KOKKOSKERNELS_BITUTILS_HPP
#include "Kokkos_Core.hpp"

#if defined(KOKKOS_COMPILER_MSVC)
#include <intrin.h>
#endif

namespace KokkosKernels {

namespace Impl {

// POP COUNT function returns the number of set bits
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned i) { return __popc(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned long i) { return __popcll(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned long long i) { return __popcll(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(int i) { return __popc(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(long i) { return __popcll(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(long long i) { return __popcll(i); }

#elif defined(__INTEL_COMPILER)
KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned i) { return _popcnt32(i); }
KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned long i) { return _popcnt64(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned long long i) { return _popcnt64(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(int i) { return _popcnt32(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(long i) { return _popcnt64(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(long long i) { return _popcnt64(i); }

#elif defined(__GNUC__) || defined(__GNUG__)
KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned i) { return __builtin_popcount(i); }
KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned long i) { return __builtin_popcountl(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned long long i) { return __builtin_popcountll(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(int i) { return __builtin_popcount(i); }
KOKKOS_FORCEINLINE_FUNCTION
int pop_count(long i) { return __builtin_popcountl(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(long long i) { return __builtin_popcountll(i); }

#elif defined(__ibmxl_vrm__)
// See
// https://www.ibm.com/support/knowledgecenter/SSGH3R_16.1.0/com.ibm.xlcpp161.aix.doc/compiler_ref/compiler_builtins.html
// link gives info about builtin names for xlclang++
KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned i) { return __builtin_popcnt4(i); }
KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned long i) { return __builtin_popcnt8(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned long long i) { return __builtin_popcnt8(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(int i) { return __builtin_popcnt4(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(long i) { return __builtin_popcnt8(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(long long i) { return __builtin_popcnt8(i); }

#elif defined(__IBMCPP__) || defined(__IBMC__)
KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned i) { return __popcnt4(i); }
KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned long i) { return __popcnt8(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned long long i) { return __popcnt8(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(int i) { return __popcnt4(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(long i) { return __popcnt8(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(long long i) { return __popcnt8(i); }

#elif defined(KOKKOS_COMPILER_MSVC)
KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned i) { return __popcnt(i); }
KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned long i) { return __popcnt(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(unsigned long long i) { return __popcnt64(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(int i) { return __popcnt(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(long i) { return __popcnt(i); }

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(long long i) { return __popcnt64(i); }

#else
#error \
    "Popcount function is not defined for this compiler. Please report this with the compiler you are using to KokkosKernels."
#endif

// least_set_bit function returns the position of right most set bit

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(unsigned i) { return __ffs(i); }

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(unsigned long i) {
#if defined(__HIP_DEVICE_COMPILE__)
  return __ffsll(static_cast<unsigned long long>(i));
#else
  return __ffsll(i);
#endif
}

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(unsigned long long i) { return __ffsll(i); }

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(int i) { return __ffs(i); }

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(long i) {
#if defined(__HIP_DEVICE_COMPILE__)
  return __ffsll(static_cast<long long>(i));
#else
  return __ffsll(i);
#endif
}

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(long long i) { return __ffsll(i); }

/*
#elif defined ( __INTEL_COMPILER )
KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit( unsigned i ){
  return _bit_scan_forward(i) + 1;
}


KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit( unsigned long long i ){
  const int llsize = sizeof(unsigned long long) * 8;
  const int intsize = sizeof(int) * 8;
  const int iteration = llsize / intsize;
  unsigned long long tmp = i;

  for (int j = 0; j < iteration; ++j){
    unsigned castint = (tmp >> (intsize * j));
    if (castint) return least_set_bit(castint) + intsize * j;
  }
  return -1;
}

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit( unsigned long i ){
  return least_set_bit( unsigned long long(i) );
}



KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(int i ){
  return _bit_scan_forward(i) + 1;
}

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit( long i ){
  return least_set_bit( unsigned long long(i) );
}


KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit( long long i ){
  return least_set_bit( unsigned long long(i) );
}
*/

#elif defined(__INTEL_COMPILER) || defined(KOKKOS_COMPILER_IBM) || defined(__GNUC__) || defined(__GNUG__)
KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(unsigned i) { return __builtin_ffs(i); }
KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(unsigned long i) { return __builtin_ffsl(i); }

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(unsigned long long i) { return __builtin_ffsll(i); }

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(int i) { return __builtin_ffs(i); }
KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(long i) { return __builtin_ffsl(i); }

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(long long i) { return __builtin_ffsll(i); }

#elif defined(KOKKOS_COMPILER_MSVC)
KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(unsigned i) { return __lzcnt(i); }
KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(unsigned long i) { return __lzcnt(i); }

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(unsigned long long i) { return __lzcnt64(i); }

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(int i) { return __lzcnt(i); }
KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(long i) { return __lzcnt(i); }

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(long long i) { return __lzcnt64(i); }

#else
#error \
    "least_set_bit function is not defined for this compiler. Please report this with the compiler you are using to KokkosKernels."
#endif

}  // namespace Impl
}  // namespace KokkosKernels

#endif

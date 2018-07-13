/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
 */

#ifndef _KOKKOSKERNELS_BITUTILS_HPP
#define _KOKKOSKERNELS_BITUTILS_HPP
#include "Kokkos_Core.hpp"

namespace KokkosKernels{

namespace Impl{

// POP COUNT function returns the number of set bits
#if defined( __CUDA_ARCH__ )
KOKKOS_FORCEINLINE_FUNCTION
int pop_count( unsigned i ){
  return __popc(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int pop_count( unsigned long i ){
  return __popcll(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int pop_count( unsigned long long i ){
  return __popcll(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int pop_count( int i ){
  return __popc(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int pop_count( long i ){
  return __popcll(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int pop_count( long long i ){
  return __popcll(i);
}

#elif defined ( __INTEL_COMPILER )
KOKKOS_FORCEINLINE_FUNCTION
int pop_count( unsigned i ){
  return _popcnt32(i);
}
KOKKOS_FORCEINLINE_FUNCTION
int pop_count( unsigned long i ){
  return _popcnt64(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int pop_count( unsigned long long i ){
  return _popcnt64(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(int i ){
  return _popcnt32(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int pop_count( long i ){
  return _popcnt64(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int pop_count( long long i ){
  return _popcnt64(i);
}

#elif defined( KOKKOS_COMPILER_IBM )
KOKKOS_FORCEINLINE_FUNCTION
int pop_count( unsigned i ){
  return __popcnt4(i);
}
KOKKOS_FORCEINLINE_FUNCTION
int pop_count( unsigned long i ){
  return __popcnt8(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int pop_count( unsigned long long i ){
  return __popcnt8(i);
}




KOKKOS_FORCEINLINE_FUNCTION
int pop_count( int i ){
  return __popcnt4(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int pop_count( long i ){
  return __popcnt8(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int pop_count( long long i ){
  return __popcnt8(i);
}

#elif defined( __GNUC__ ) || defined( __GNUG__ )
KOKKOS_FORCEINLINE_FUNCTION
int pop_count( unsigned i ){
  return __builtin_popcount(i);
}
KOKKOS_FORCEINLINE_FUNCTION
int pop_count( unsigned long i ){
  return __builtin_popcountl(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int pop_count( unsigned long long i ){
  return __builtin_popcountll(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int pop_count( int i ){
  return __builtin_popcount(i);
}
KOKKOS_FORCEINLINE_FUNCTION
int pop_count(  long i ){
  return __builtin_popcountl(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int pop_count(  long long i ){
  return __builtin_popcountll(i);
}

#else
  #error "Popcount function is not defined for this compiler. Please report this with the compiler you are using to KokkosKernels."
#endif


// least_set_bit function returns the position of right most set bit

#if defined( __CUDA_ARCH__ )
KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit( unsigned i ){
  return __ffs(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit( unsigned long i ){
  return __ffsll(i);
}


KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit( unsigned long long i ){
  return __ffsll(i);
}



KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit( int i ){
  return __ffs(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit( long i ){
  return __ffsll(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit( long long i ){
  return __ffsll(i);
}


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

#elif defined ( __INTEL_COMPILER ) || defined( KOKKOS_COMPILER_IBM ) || defined( __GNUC__ ) || defined( __GNUG__ )
KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit( unsigned i ){
  return __builtin_ffs(i);
}
KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit( unsigned long i ){
  return __builtin_ffsl(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit( unsigned long long i ){
  return __builtin_ffsll(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit( int i ){
  return __builtin_ffs(i);
}
KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(  long i ){
  return __builtin_ffsl(i);
}

KOKKOS_FORCEINLINE_FUNCTION
int least_set_bit(  long long i ){
  return __builtin_ffsll(i);
}

#else
  #error "least_set_bit function is not defined for this compiler. Please report this with the compiler you are using to KokkosKernels."
#endif


}
}

#endif

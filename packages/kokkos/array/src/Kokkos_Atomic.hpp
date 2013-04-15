/*
//@HEADER
// ************************************************************************
//
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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


#ifndef KOKKOSARRAY_ATOMIC_H_
#define KOKKOSARRAY_ATOMIC_H_


//This header file defines prototypes for the atomic functions: exchange, compare_exchange, add
//Supported types are: signed and unsigned 4 and 8 byte integers, float, and double
//They are implemented through GCC compatible intrinsics, OpenMP directives and native CUDA intrinsics.

#include <stdint.h>
#include "KokkosArray_Macros.hpp"
#include "KokkosArray_Host.hpp"

namespace KokkosArray {

#if defined(GNU_ATOMICS_GCC) || defined(GNU_ATOMICS_INTEL) || defined(OMP31_ATOMICS)
#define EXTERNAL_ATOMICS_CHOICE
#endif

#ifndef EXTERNAL_ATOMICS_CHOICE
#ifndef __CUDA_ARCH__
  #undef EXTERNAL_ATOMICS_CHOICE
  #if !defined(_OPENMP) || _OPENMP < 201107
    #if defined(__GNUC__) || defined(__GNUG__)
      #define GNU_ATOMICS_GCC
    #else
      #define GNU_ATOMICS_INTEL
    #endif
  #else
    #define OMP31_ATOMICS
  #endif
#else
  #define CUDA_ATOMICS
#endif
#endif


//replaces the value at dest with val, and returns the prior value.
//template<typename T>
//KOKKOSARRAY_INLINE_FUNCTION T atomic_exchange(volatile T* dest, T val);

#include "impl/Kokkos_Atomic_Exchange.hpp"



#ifndef EXTERNAL_ATOMICS_CHOICE
#undef GNU_ATOMICS_GCC
#undef GNU_ATOMICS_INTEL
#undef OMP31_ATOMICS
#undef CUDA_ATOMICS
#ifndef __CUDA_ARCH__
  #if !defined(_OPENMP) || _OPENMP < 201107
    #if defined(__GNUC__) || defined(__GNUG__)
      #define GNU_ATOMICS_GCC
    #else
      #define GNU_ATOMICS_INTEL
    #endif
  #else
    #define OMP31_ATOMICS
  #endif
#else
  #define CUDA_ATOMICS
#endif
#endif

const char* atomic_query_version() {
#ifdef GNU_ATOMICS_GCC
	return "GNU_ATOMICS_GCC";
#endif
#ifdef GNU_ATOMICS_INTEL
	return "GNU_ATOMICS_INTEL";
#endif
#ifdef OMP31_ATOMICS
	return "OMP31_ATOMICS";
#endif
#ifdef CUDA_ATOMICS
	return "CUDA_ATOMICS";
#endif
}

//template<class T>
//bool atomic_compare_exchange_strong(volatile T* dest, T compare, T val);
#include "impl/Kokkos_Atomic_Compare_Exchange_Strong.hpp"

#ifndef EXTERNAL_ATOMICS_CHOICE
#undef GNU_ATOMICS_GCC
#undef GNU_ATOMICS_INTEL
#undef OMP31_ATOMICS
#undef CUDA_ATOMICS
#ifndef __CUDA_ARCH__
  #if !defined(_OPENMP) || _OPENMP < 201107
    #if defined(__GNUC__) || defined(__GNUG__)
      #define GNU_ATOMICS_GCC
    #else
      #define GNU_ATOMICS_INTEL
    #endif
  #else
    #define OMP31_ATOMICS
  #endif
#else
  #define CUDA_ATOMICS
#endif
#endif

//template<class T>
//T atomic_fetch_add(volatile T* dest, T val);
#include "impl/Kokkos_Atomic_Fetch_Add.hpp"

}

#endif /* KOKKOSARRAY_ATOMIC_H_ */

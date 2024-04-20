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

#ifndef STK_SIMD_INDEX_H
#define STK_SIMD_INDEX_H

#include <stk_simd/Simd.hpp>

namespace stk { namespace simd {

struct Index {
  KOKKOS_INLINE_FUNCTION explicit Index(int i) : index(i) {}
  friend KOKKOS_INLINE_FUNCTION int int_index(const Index&);
 private:
  int index;
};

KOKKOS_INLINE_FUNCTION int int_index(const Index& i) {
  return i.index;
}

KOKKOS_INLINE_FUNCTION int int_index(const int& i) {
  return i;
}

template <typename T>
struct IndexTraits {
  typedef simd::Double double_type;
  typedef simd::Float float_type;
};

template <>
struct IndexTraits<int> {
  typedef double double_type;
  typedef float float_type;
};

#if defined(STK_ENABLE_GPU) || defined(USE_STK_SIMD_NONE)
typedef int DeviceIndex;

template <typename T>
struct DeviceTraits {
  typedef typename stk::Traits<T>::base_type simd_type;
};
#else
typedef Index DeviceIndex;

template <typename T>
struct DeviceTraits {
  typedef typename stk::Traits<T>::simd_type simd_type;
};
#endif

typedef typename IndexTraits<DeviceIndex>::double_type DeviceDouble;
typedef typename IndexTraits<DeviceIndex>::float_type  DeviceFloat;
}}



#endif

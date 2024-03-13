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

#ifndef STK_UTIL_UTIL_STRIDEDARRAY_HPP
#define STK_UTIL_UTIL_STRIDEDARRAY_HPP

#include <stk_util/stk_config.h>
#include <stk_util/util/PairIter.hpp>
#include <Kokkos_Core.hpp>

namespace stk
{
namespace util
{
constexpr unsigned defaultStride = 512;

template <typename T>
class StridedArray
{
public:
  KOKKOS_FUNCTION 
  StridedArray()
  : dataPointer(nullptr),
    length(0)
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    , stride(defaultStride)
#endif
  {
  }

  KOKKOS_FUNCTION
  StridedArray(T* e,
               unsigned n,
               int stride_in=defaultStride)
  : dataPointer(e),
    length(n)
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    , stride(stride_in)
#endif
  {
  }

  KOKKOS_FUNCTION
  StridedArray(PairIter<T*> data
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
               , int stride_in=defaultStride
#endif
              )
  : dataPointer(data.begin()),
    length(data.size())
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    , stride(stride_in)
#endif
  {
  }

  KOKKOS_FUNCTION
  T operator[](unsigned i) const
  {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    return dataPointer[stride*i];
#else
    return dataPointer[i];
#endif
  }

  KOKKOS_FUNCTION
  T* data() { return dataPointer; }
  KOKKOS_FUNCTION
  const T* data() const { return dataPointer; }

  KOKKOS_FUNCTION
  unsigned size() const
  { 
    return length;
  }

private:
  T* dataPointer;
  unsigned length;
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  int stride;
#endif
};

} //namespace util
} //namespace stk

#endif

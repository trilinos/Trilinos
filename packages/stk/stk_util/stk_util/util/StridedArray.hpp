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
#include "Kokkos_Macros.hpp"
#include <iterator>
#include <cstddef>

namespace stk
{
namespace util
{
constexpr unsigned defaultStride = 1;

template <typename T>
class StridedArray
{
public:
  KOKKOS_INLINE_FUNCTION 
  StridedArray()
  : dataPointer(nullptr),
    length(0)
#ifdef STK_ENABLE_GPU
    , stride(defaultStride)
#endif
  {
  }

  KOKKOS_INLINE_FUNCTION
  StridedArray(T* e,
               unsigned n,
               [[maybe_unused]] int stride_in=defaultStride)
  : dataPointer(e),
    length(n)
#ifdef STK_ENABLE_GPU
    , stride(stride_in)
#endif
  {
  }

  KOKKOS_INLINE_FUNCTION
  StridedArray(PairIter<T*> data,
               [[maybe_unused]] int stride_in=defaultStride)
  : dataPointer(data.begin()),
    length(data.size())
#ifdef STK_ENABLE_GPU
    , stride(stride_in)
#endif
  {
  }

  KOKKOS_INLINE_FUNCTION
  T operator[](unsigned i) const
  {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    return dataPointer[stride*i];
#else
    return dataPointer[i];
#endif
  }

  KOKKOS_INLINE_FUNCTION
  T* data() { return dataPointer; }
  KOKKOS_INLINE_FUNCTION
  const T* data() const { return dataPointer; }

  KOKKOS_INLINE_FUNCTION
  unsigned size() const
  { 
    return length;
  }

  KOKKOS_INLINE_FUNCTION
  bool empty() const
  {
    return size() == 0;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator==(const StridedArray<T>& rhs) const
  {
    if (this->size() != rhs.size()) {
      return false;
    }

    for(unsigned i=0; i<size(); ++i) {
      if ((*this)[i] != rhs[i]) {
        return false;
      }
    }

    return true;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator!=(const StridedArray<T>& rhs) const
  {
    return !(*this == rhs);
  }

  // Begin and end methods for range-based for loops
  class Iterator {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;

    KOKKOS_INLINE_FUNCTION
    Iterator(T* ptr,
             [[maybe_unused]] int stride_in=defaultStride)
    : m_ptr(ptr)
#ifdef STK_ENABLE_GPU
      , m_stride(stride_in)
#endif
    {}

    KOKKOS_INLINE_FUNCTION
    const T& operator*() const { return *m_ptr; }

    KOKKOS_INLINE_FUNCTION
    Iterator& operator++() {
#ifdef STK_ENABLE_GPU
      m_ptr += m_stride;
#else
      ++m_ptr;
#endif
      return *this;
    }

    KOKKOS_INLINE_FUNCTION
    difference_type operator-(const Iterator& other) const { return m_ptr - other.m_ptr; }

    KOKKOS_INLINE_FUNCTION
    bool operator!=(const Iterator& other) const { return m_ptr != other.m_ptr; }
    KOKKOS_INLINE_FUNCTION
    bool operator==(const Iterator& other) const { return m_ptr == other.m_ptr; }

  private:
    T* m_ptr;
#ifdef STK_ENABLE_GPU
    int m_stride;
#endif
  };

  KOKKOS_INLINE_FUNCTION
  Iterator begin() const {
#ifdef STK_ENABLE_GPU
    return Iterator(dataPointer, stride);
#else
    return Iterator(dataPointer);
#endif
  }

  KOKKOS_INLINE_FUNCTION
  Iterator end() const {
#ifdef STK_ENABLE_GPU
    return Iterator(dataPointer + stride * length, stride);
#else
    return Iterator(dataPointer + length, 1);
#endif
  }

private:
  T* dataPointer;
  unsigned length;
#ifdef STK_ENABLE_GPU
  int stride;
#endif
};

} //namespace util
} //namespace stk

#endif

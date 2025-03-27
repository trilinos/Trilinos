/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

// IWYU pragma: private

#ifndef STK_INCLUDE_ONLY_STK_SIMD_HEADER
static_assert(false, "Do not include simd impl files directly. Only include stk_simd/Simd.hpp");
#endif

#pragma once

#include "simd_common.hpp"

#ifdef __CUDACC__
#define SIMD_CUDA_ALWAYS_INLINE __forceinline__
#endif

#ifdef __CUDACC__
#define SIMD_HOST_DEVICE __host__ __device__
#else
#define SIMD_HOST_DEVICE
#endif

#ifdef __CUDACC__
#define SIMD_DEVICE __device__
#else
#define SIMD_DEVICE
#endif

#ifdef __CUDACC__

namespace SIMD_NAMESPACE {

namespace simd_abi {

template <int N>
class cuda_warp {
  static_assert(N <= 32, "CUDA warps can't be more than 32 threads");
 public:
  SIMD_HOST_DEVICE static unsigned mask() {
    return (unsigned(1) << N) - unsigned(1);
  }
};

}

template <class T, int N>
class simd_storage<T, simd_abi::cuda_warp<N>> {
  T m_value[simd<T, simd_abi::cuda_warp<N>>::size()];
 public:
  using value_type = T;
  using abi_type = simd_abi::cuda_warp<N>;
  using simd_type = simd<T, abi_type>;
  SIMD_ALWAYS_INLINE inline simd_storage() = default;
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline static constexpr
  int size() { return simd<T, abi_type>::size(); }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline
  simd_storage(simd<T, abi_type> const& value) {
    value.copy_to(m_value, element_aligned_tag());
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE explicit inline
  simd_storage(T value)
    :simd_storage(simd<T, abi_type>(value))
  {}
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline
  simd_storage& operator=(simd<T, abi_type> const& value) {
    value.copy_to(m_value, element_aligned_tag());
    return *this;
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE
  T const* data() const { return m_value; }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE
  T* data() { return m_value; }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE
  T const& operator[](int i) const { return m_value[i]; }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE
  T& operator[](int i) { return m_value[i]; }
};

template <class T, int N>
class simd_mask<T, simd_abi::cuda_warp<N>> {
  bool m_value;
 public:
  using value_type = bool;
  using abi_type = simd_abi::cuda_warp<N>;
  using simd_type = simd<T, abi_type>;
  SIMD_CUDA_ALWAYS_INLINE simd_mask() = default;
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE static constexpr
  int size() { return N; }
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE
  simd_mask(bool value)
    :m_value(value)
  {}
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE constexpr
  bool get() const {
    return m_value;
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE
  simd_mask operator||(simd_mask const& other) const {
    return m_value || other.m_value;
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE
  simd_mask operator&&(simd_mask const& other) const {
    return m_value && other.m_value;
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE
  simd_mask operator!() const {
    return !m_value;
  }
};

template <class T, int N>
SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE
bool all_of(simd_mask<T, simd_abi::cuda_warp<N>> const& a) {
  return bool(__all_sync(simd_abi::cuda_warp<N>::mask(), int(a.get())));
}

template <class T, int N>
SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE
bool any_of(simd_mask<T, simd_abi::cuda_warp<N>> const& a) {
  return bool(__any_sync(simd_abi::cuda_warp<N>::mask(), int(a.get())));
}

template <class T, int N>
class simd<T, simd_abi::cuda_warp<N>> {
  T m_value;
 public:
  using value_type = T;
  using abi_type = simd_abi::cuda_warp<N>;
  using mask_type = simd_mask<T, abi_type>;
  using storage_type = simd_storage<T, abi_type>;
  SIMD_CUDA_ALWAYS_INLINE simd() = default;
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE static constexpr int size() { return N; }
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE simd(T value)
    :m_value(value)
  {}
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE
  simd(storage_type const& value) {
    copy_from(value.data(), element_aligned_tag());
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE
  simd& operator=(storage_type const& value) {
    copy_from(value.data(), element_aligned_tag());
    return *this;
  }
  template <class Flags>
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE simd(T const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE simd operator*(simd const& other) const {
    return simd(m_value * other.m_value);
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE simd operator/(simd const& other) const {
    return simd(m_value / other.m_value);
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE simd operator+(simd const& other) const {
    return simd(m_value + other.m_value);
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE simd operator-(simd const& other) const {
    return simd(m_value - other.m_value);
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE simd operator-() const {
    return simd(-m_value);
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE void copy_from(T const* ptr, element_aligned_tag) {
#ifdef __CUDA_ARCH__
    m_value = ptr[threadIdx.x];
#endif
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE void copy_to(T* ptr, element_aligned_tag) const {
#ifdef __CUDA_ARCH__
    ptr[threadIdx.x] = m_value;
#endif
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE T get() const {
    return m_value;
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE
  mask_type operator<(simd const& other) const {
    return mask_type(m_value < other.m_value);
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE
  mask_type operator==(simd const& other) const {
    return mask_type(m_value == other.m_value);
  }
};

template <class T, int N>
SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE simd<T, simd_abi::cuda_warp<N>> abs(simd<T, simd_abi::cuda_warp<N>> const& a) {
  return simd<T, simd_abi::cuda_warp<N>>(std::abs(a.get()));
}

template <class T, int N>
SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE simd<T, simd_abi::cuda_warp<N>> sqrt(simd<T, simd_abi::cuda_warp<N>> const& a) {
  return simd<T, simd_abi::cuda_warp<N>>(std::sqrt(a.get()));
}

template <class T, int N>
SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE simd<T, simd_abi::cuda_warp<N>> cbrt(simd<T, simd_abi::cuda_warp<N>> const& a) {
  return simd<T, simd_abi::cuda_warp<N>>(std::cbrt(a.get()));
}

template <class T, int N>
SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE simd<T, simd_abi::cuda_warp<N>> exp(simd<T, simd_abi::cuda_warp<N>> const& a) {
  return simd<T, simd_abi::cuda_warp<N>>(std::exp(a.get()));
}

template <class T, int N>
SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE simd<T, simd_abi::cuda_warp<N>> fma(
    simd<T, simd_abi::cuda_warp<N>> const& a,
    simd<T, simd_abi::cuda_warp<N>> const& b,
    simd<T, simd_abi::cuda_warp<N>> const& c) {
  return simd<T, simd_abi::cuda_warp<N>>((a.get() * b.get()) + c.get());
}

template <class T, int N>
SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE simd<T, simd_abi::cuda_warp<N>> max(
    simd<T, simd_abi::cuda_warp<N>> const& a, simd<T, simd_abi::cuda_warp<N>> const& b) {
  return simd<T, simd_abi::cuda_warp<N>>((a.get() < b.get()) ? b.get() : a.get());
}

template <class T, int N>
SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE simd<T, simd_abi::cuda_warp<N>> min(
    simd<T, simd_abi::cuda_warp<N>> const& a, simd<T, simd_abi::cuda_warp<N>> const& b) {
  return simd<T, simd_abi::cuda_warp<N>>((b.get() < a.get()) ? b.get() : a.get());
}

template <class T, int N>
SIMD_CUDA_ALWAYS_INLINE SIMD_HOST_DEVICE simd<T, simd_abi::cuda_warp<N>> choose(
    simd_mask<T, simd_abi::cuda_warp<N>> const& a,
    simd<T, simd_abi::cuda_warp<N>> const& b,
    simd<T, simd_abi::cuda_warp<N>> const& c) {
  return simd<T, simd_abi::cuda_warp<N>>(a.get() ? b.get() : c.get());
}

}

#endif

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

namespace SIMD_NAMESPACE {

namespace simd_abi {

template <int N>
class pack;

}

template <int N>
class simd_mask<float, simd_abi::pack<N>> {
  int m_value[N];
 public:
  using value_type = bool;
  using simd_type = simd<float, simd_abi::pack<N>>;
  using abi_type = simd_abi::pack<N>;
  SIMD_ALWAYS_INLINE inline simd_mask() = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return N; }
  SIMD_ALWAYS_INLINE inline simd_mask(bool value) {
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) m_value[i] = value;
  }
  SIMD_ALWAYS_INLINE inline constexpr bool operator[](int i) const { return m_value[i]; }
  SIMD_ALWAYS_INLINE inline int& operator[](int i) { return m_value[i]; }
  SIMD_ALWAYS_INLINE inline simd_mask operator||(simd_mask const& other) const {
    simd_mask result;
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) result.m_value[i] = m_value[i] || other.m_value[i];
    return result;
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator&&(simd_mask const& other) const {
    simd_mask result;
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) result.m_value[i] = m_value[i] && other.m_value[i];
    return result;
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator!() const {
    simd_mask result;
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) result.m_value[i] = !m_value[i];
    return result;
  }
};

template <int N>
class simd_mask<double, simd_abi::pack<N>> {
  std::int64_t m_value[N];
 public:
  using value_type = bool;
  using simd_type = simd<double, simd_abi::pack<N>>;
  using abi_type = simd_abi::pack<N>;
  SIMD_ALWAYS_INLINE inline simd_mask() = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return N; }
  SIMD_ALWAYS_INLINE inline simd_mask(bool value) {
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) m_value[i] = value;
  }
  SIMD_ALWAYS_INLINE inline constexpr bool operator[](int i) const { return m_value[i]; }
  SIMD_ALWAYS_INLINE inline std::int64_t& operator[](int i) { return m_value[i]; }
  SIMD_ALWAYS_INLINE inline simd_mask operator||(simd_mask const& other) const {
    simd_mask result;
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) result.m_value[i] = m_value[i] || other.m_value[i];
    return result;
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator&&(simd_mask const& other) const {
    simd_mask result;
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) result.m_value[i] = m_value[i] && other.m_value[i];
    return result;
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator!() const {
    simd_mask result;
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) result.m_value[i] = !m_value[i];
    return result;
  }
};

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline bool all_of(simd_mask<T, simd_abi::pack<N>> const& a) {
  bool result = true;
  SIMD_PRAGMA for (int i = 0; i < a.size(); ++i) result = result && a[i];
  return result;
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline bool any_of(simd_mask<T, simd_abi::pack<N>> const& a) {
  bool result = false;
  SIMD_PRAGMA for (int i = 0; i < a.size(); ++i) result = result || a[i];
  return result;
}

template <class T, int N>
class simd<T, simd_abi::pack<N>> {
  T m_value[N];
 public:
  using value_type = T;
  using abi_type = simd_abi::pack<N>;
  using mask_type = simd_mask<T, abi_type>;
  using storage_type = simd_storage<T, abi_type>;
  SIMD_ALWAYS_INLINE inline simd() = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return N; }
  SIMD_ALWAYS_INLINE inline simd(T value)
  {
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) m_value[i] = value;
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline
  simd(storage_type const& value) {
    copy_from(value.data(), element_aligned_tag());
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline
  simd& operator=(storage_type const& value) {
    copy_from(value.data(), element_aligned_tag());
    return *this;
  }
  template <class Flags>
  SIMD_ALWAYS_INLINE simd(T const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_ALWAYS_INLINE simd operator*(simd const& other) const {
    simd result;
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) result[i] = m_value[i] * other.m_value[i];
    return result;
  }
  SIMD_ALWAYS_INLINE simd operator/(simd const& other) const {
    simd result;
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) result[i] = m_value[i] / other.m_value[i];
    return result;
  }
  SIMD_ALWAYS_INLINE simd operator+(simd const& other) const {
    simd result;
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) result[i] = m_value[i] + other.m_value[i];
    return result;
  }
  SIMD_ALWAYS_INLINE simd operator-(simd const& other) const {
    simd result;
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) result[i] = m_value[i] - other.m_value[i];
    return result;
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd operator-() const {
    simd result;
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) result[i] = -m_value[i];
    return result;
  }
  SIMD_ALWAYS_INLINE void copy_from(T const* ptr, element_aligned_tag) {
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) m_value[i] = ptr[i];
  }
  SIMD_ALWAYS_INLINE void copy_to(T* ptr, element_aligned_tag) const {
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) ptr[i] = m_value[i];
  }
  SIMD_ALWAYS_INLINE constexpr T operator[](int i) const { return m_value[i]; }
  SIMD_ALWAYS_INLINE T& operator[](int i) { return m_value[i]; }
  SIMD_ALWAYS_INLINE simd_mask<T, simd_abi::pack<N>> operator<(simd const& other) const {
    simd_mask<T, simd_abi::pack<N>> result;
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) result[i] = m_value[i] < other.m_value[i];
    return result;
  }
  SIMD_ALWAYS_INLINE simd_mask<T, simd_abi::pack<N>> operator==(simd const& other) const {
    simd_mask<T, simd_abi::pack<N>> result;
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) result[i] = m_value[i] == other.m_value[i];
    return result;
  }
};

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::pack<N>> abs(simd<T, simd_abi::pack<N>> const& a) {
  simd<T, simd_abi::pack<N>> result;
  using std::abs;
  SIMD_PRAGMA for (int i = 0; i < a.size(); ++i) result[i] = abs(a[i]);
  return result;
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::pack<N>> sqrt(simd<T, simd_abi::pack<N>> const& a) {
  simd<T, simd_abi::pack<N>> result;
  using std::sqrt;
  SIMD_PRAGMA for (int i = 0; i < a.size(); ++i) result[i] = sqrt(a[i]);
  return result;
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::pack<N>> cbrt(simd<T, simd_abi::pack<N>> const& a) {
  simd<T, simd_abi::pack<N>> result;
  using std::cbrt;
  SIMD_PRAGMA for (int i = 0; i < a.size(); ++i) result[i] = cbrt(a[i]);
  return result;
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::pack<N>> exp(simd<T, simd_abi::pack<N>> const& a) {
  simd<T, simd_abi::pack<N>> result;
  using std::exp;
  SIMD_PRAGMA for (int i = 0; i < a.size(); ++i) result[i] = exp(a[i]);
  return result;
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::pack<N>> fma(
    simd<T, simd_abi::pack<N>> const& a,
    simd<T, simd_abi::pack<N>> const& b,
    simd<T, simd_abi::pack<N>> const& c) {
  simd<T, simd_abi::pack<N>> result;
  SIMD_PRAGMA for (int i = 0; i < a.size(); ++i) result[i] = (a[i] * b[i]) + c[i];
  return result;
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::pack<N>> max(
    simd<T, simd_abi::pack<N>> const& a, simd<T, simd_abi::pack<N>> const& b) {
  simd<T, simd_abi::pack<N>> result;
  SIMD_PRAGMA
  for (int i = 0; i < a.size(); ++i) {
    result[i] = choose((a[i] < b[i]), b[i], a[i]);
  }
  return result;
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::pack<N>> min(
    simd<T, simd_abi::pack<N>> const& a, simd<T, simd_abi::pack<N>> const& b) {
  simd<T, simd_abi::pack<N>> result;
  SIMD_PRAGMA
  for (int i = 0; i < a.size(); ++i) {
    result[i] = choose((b[i] < a[i]), b[i], a[i]);
  }
  return result;
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::pack<N>> choose(
    simd_mask<T, simd_abi::pack<N>> const& a, simd<T, simd_abi::pack<N>> const& b, simd<T, simd_abi::pack<N>> const& c) {
  simd<T, simd_abi::pack<N>> result;
  SIMD_PRAGMA for (int i = 0; i < a.size(); ++i) result[i] = a[i] ? b[i] : c[i];
  return result;
}

}

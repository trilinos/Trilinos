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

#if defined(__clang__)

namespace SIMD_NAMESPACE {

namespace simd_abi {

template <int N>
class vector_size {};

}

template <int N>
class simd_mask<float, simd_abi::vector_size<N>> {
  typedef int native_type __attribute__((vector_size(N)));
  native_type m_value;
 public:
  using value_type = bool;
  using simd_type = simd<float, simd_abi::vector_size<N>>;
  using abi_type = simd_abi::vector_size<N>;
  SIMD_ALWAYS_INLINE inline simd_mask() = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return N / sizeof(int); }
  SIMD_ALWAYS_INLINE inline simd_mask(bool value)
    :m_value(-int(value))
  {}
  SIMD_ALWAYS_INLINE inline simd_mask(native_type value)
    :m_value(value)
  {}
  SIMD_ALWAYS_INLINE inline int operator[](int i) { return m_value[i]; }
  SIMD_ALWAYS_INLINE inline native_type const& get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask operator||(simd_mask const& other) const {
    return simd_mask(m_value || other.m_value);
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator&&(simd_mask const& other) const {
    return simd_mask(m_value && other.m_value);
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator!() const {
    return simd_mask(!m_value);
  }
};

template <int N>
class simd_mask<double, simd_abi::vector_size<N>> {
  typedef long long native_type __attribute__((vector_size(N)));
  native_type m_value;
 public:
  using value_type = bool;
  using simd_type = simd<double, simd_abi::vector_size<N>>;
  using abi_type = simd_abi::vector_size<N>;
  SIMD_ALWAYS_INLINE inline simd_mask() = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return N / sizeof(long long); }
  SIMD_ALWAYS_INLINE inline simd_mask(bool value);
  SIMD_ALWAYS_INLINE inline simd_mask(native_type value)
    :m_value(value)
  {}
  SIMD_ALWAYS_INLINE inline long long operator[](int i) { return m_value[i]; }
  SIMD_ALWAYS_INLINE inline native_type const& get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask operator||(simd_mask const& other) const {
    return simd_mask(m_value || other.m_value);
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator&&(simd_mask const& other) const {
    return simd_mask(m_value && other.m_value);
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator!() const {
    return simd_mask(!m_value);
  }
};

template <>
SIMD_ALWAYS_INLINE inline simd_mask<float, simd_abi::vector_size<32>>::simd_mask(bool value)
{
  m_value = {-int(value), -int(value), -int(value), -int(value),
             -int(value), -int(value), -int(value), -int(value)};
}

template <>
SIMD_ALWAYS_INLINE inline simd_mask<double, simd_abi::vector_size<32>>::simd_mask(bool value)
{
  m_value = {-(long long)(value), -(long long)(value), -(long long)(value), -(long long)(value)};
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline bool all_of(simd_mask<T, simd_abi::vector_size<N>> const& a) {
  bool result = true;
  SIMD_PRAGMA for (int i = 0; i < a.size(); ++i) result = result && a.get()[i];
  return result;
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline bool any_of(simd_mask<T, simd_abi::vector_size<N>> const& a) {
  bool result = false;
  SIMD_PRAGMA for (int i = 0; i < a.size(); ++i) result = result || a.get()[i];
  return result;
}

template <class T, int N>
class simd<T, simd_abi::vector_size<N>> {
  typedef T native_type __attribute__((vector_size(N)));
  native_type m_value;
 public:
  using value_type = T;
  using abi_type = simd_abi::vector_size<N>;
  using mask_type = simd_mask<T, abi_type>;
  using storage_type = simd_storage<T, abi_type>;
  SIMD_ALWAYS_INLINE inline simd() = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return N / sizeof(T); }
  SIMD_ALWAYS_INLINE inline simd(T value);
  SIMD_ALWAYS_INLINE inline simd(native_type value):m_value(value) {}
  SIMD_ALWAYS_INLINE inline
  simd(storage_type const& value) {
    copy_from(value.data(), element_aligned_tag());
  }
  SIMD_ALWAYS_INLINE inline
  simd& operator=(storage_type const& value) {
    copy_from(value.data(), element_aligned_tag());
    return *this;
  }
  template <class Flags>
  SIMD_ALWAYS_INLINE simd(T const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_ALWAYS_INLINE simd operator*(simd const& other) const {
    return simd(m_value * other.m_value);
  }
  SIMD_ALWAYS_INLINE simd operator/(simd const& other) const {
    return simd(m_value / other.m_value);
  }
  SIMD_ALWAYS_INLINE simd operator+(simd const& other) const {
    return simd(m_value + other.m_value);
  }
  SIMD_ALWAYS_INLINE simd operator-(simd const& other) const {
    return simd(m_value - other.m_value);
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd operator-() const {
    return simd(-m_value);
  }
  SIMD_ALWAYS_INLINE void copy_from(T const* ptr, element_aligned_tag) {
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) m_value[i] = ptr[i];
  }
  SIMD_ALWAYS_INLINE void copy_to(T* ptr, element_aligned_tag) const {
    SIMD_PRAGMA for (int i = 0; i < size(); ++i) ptr[i] = m_value[i];
  }
  SIMD_ALWAYS_INLINE constexpr T operator[](int i) const { return m_value[i]; }
  SIMD_ALWAYS_INLINE native_type const& get() const { return m_value; }
  SIMD_ALWAYS_INLINE native_type& get() { return m_value; }
  SIMD_ALWAYS_INLINE simd_mask<T, simd_abi::vector_size<N>> operator<(simd const& other) const {
    return simd_mask<T, simd_abi::vector_size<N>>(m_value < other.m_value);
  }
  SIMD_ALWAYS_INLINE simd_mask<T, simd_abi::vector_size<N>> operator==(simd const& other) const {
    return simd_mask<T, simd_abi::vector_size<N>>(m_value == other.m_value);
  }
};

template <>
SIMD_ALWAYS_INLINE inline simd<float, simd_abi::vector_size<32>>::simd(float value) {
  m_value = {value, value, value, value,
             value, value, value, value};
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::vector_size<N>> abs(simd<T, simd_abi::vector_size<N>> const& a) {
  simd<T, simd_abi::vector_size<N>> result;
  using std::sqrt;
  SIMD_PRAGMA for (int i = 0; i < a.size(); ++i) result.get()[i] = abs(a[i]);
  return result;
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::vector_size<N>> sqrt(simd<T, simd_abi::vector_size<N>> const& a) {
  simd<T, simd_abi::vector_size<N>> result;
  using std::sqrt;
  SIMD_PRAGMA for (int i = 0; i < a.size(); ++i) result.get()[i] = sqrt(a[i]);
  return result;
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::vector_size<N>> cbrt(simd<T, simd_abi::vector_size<N>> const& a) {
  simd<T, simd_abi::vector_size<N>> result;
  using std::cbrt;
  SIMD_PRAGMA for (int i = 0; i < a.size(); ++i) result.get()[i] = cbrt(a[i]);
  return result;
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::vector_size<N>> exp(simd<T, simd_abi::vector_size<N>> const& a) {
  simd<T, simd_abi::vector_size<N>> result;
  using std::exp;
  SIMD_PRAGMA for (int i = 0; i < a.size(); ++i) result.get()[i] = exp(a[i]);
  return result;
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::vector_size<N>> fma(
    simd<T, simd_abi::vector_size<N>> const& a,
    simd<T, simd_abi::vector_size<N>> const& b,
    simd<T, simd_abi::vector_size<N>> const& c) {
  simd<T, simd_abi::vector_size<N>> result;
  SIMD_PRAGMA for (int i = 0; i < a.size(); ++i) result.get()[i] = (a[i] * b[i]) + c[i];
  return result;
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::vector_size<N>> choose(
    simd_mask<T, simd_abi::vector_size<N>> const& a,
    simd<T, simd_abi::vector_size<N>> const& b,
    simd<T, simd_abi::vector_size<N>> const& c) {
  simd<T, simd_abi::vector_size<N>> result;
  SIMD_PRAGMA for (int i = 0; i < a.size(); ++i) result.get()[i] = a.get()[i] ? b.get()[i] : c.get()[i];
  return result;
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::vector_size<N>> max(
    simd<T, simd_abi::vector_size<N>> const& a,
    simd<T, simd_abi::vector_size<N>> const& b) {
  return choose(b < a, a, b);
}

template <class T, int N>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::vector_size<N>> min(
    simd<T, simd_abi::vector_size<N>> const& a,
    simd<T, simd_abi::vector_size<N>> const& b) {
  return choose(a < b, a, b);
}

}

#endif

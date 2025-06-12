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

class scalar {};

}

template <class T>
class simd_mask<T, simd_abi::scalar> {
  bool m_value;
 public:
  using value_type = bool;
  using simd_type = simd<T, simd_abi::scalar>;
  using abi_type = simd_abi::scalar;
  SIMD_ALWAYS_INLINE inline simd_mask() = default;
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE static constexpr int size() { return 1; }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd_mask(bool value)
    :m_value(value)
  {}
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline constexpr bool get() const { return m_value; }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd_mask operator||(simd_mask const& other) const {
    return m_value || other.m_value;
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd_mask operator&&(simd_mask const& other) const {
    return m_value && other.m_value;
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd_mask operator!() const {
    return !m_value;
  }
};

template <class T>
class simd_storage<T, simd_abi::scalar> {
  using Abi = simd_abi::scalar;
  T m_value;
 public:
  using value_type = T;
  using simd_type = simd<T, Abi>;
  SIMD_ALWAYS_INLINE inline simd_storage() = default;
  SIMD_ALWAYS_INLINE inline static constexpr
  int size() { return simd<T, Abi>::size(); }
  SIMD_ALWAYS_INLINE explicit SIMD_HOST_DEVICE
  simd_storage(simd<T, Abi> const& value) SIMD_HOST_DEVICE {
    value.copy_to(&m_value, element_aligned_tag());
  }
  SIMD_ALWAYS_INLINE explicit SIMD_HOST_DEVICE
  simd_storage(T value)
    :simd_storage(simd<T, Abi>(value))
  {}
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE
  simd_storage& operator=(simd<T, Abi> const& value) {
    value.copy_to(&m_value, element_aligned_tag());
    return *this;
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE
  T const* data() const { return &m_value; }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE
  T* data() { return &m_value; }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE
  T const& operator[](int) const { return m_value; }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE
  T& operator[](int) { return m_value; }
};

template <class T>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline
bool all_of(simd_mask<T, simd_abi::scalar> const& a) { return a.get(); }

template <class T>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline
bool any_of(simd_mask<T, simd_abi::scalar> const& a) { return a.get(); }

template <class T>
class simd<T, simd_abi::scalar> {
  T m_value;
 public:
  using value_type = T;
  using abi_type = simd_abi::scalar;
  using mask_type = simd_mask<T, abi_type>;
  using storage_type = simd_storage<T, abi_type>;
  SIMD_ALWAYS_INLINE inline simd() = default;
  SIMD_ALWAYS_INLINE inline simd(simd const&) = default;
  SIMD_ALWAYS_INLINE inline simd(simd&&) = default;
  SIMD_ALWAYS_INLINE inline simd& operator=(simd const&) = default;
  SIMD_ALWAYS_INLINE inline simd& operator=(simd&&) = default;
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE static constexpr int size() { return 1; }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd(T value)
    :m_value(value)
  {}
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline
  simd(storage_type const& value) {
    copy_from(value.data(), element_aligned_tag());
  }
#ifdef STK_VOLATILE_SIMD
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline
  simd(simd const volatile& value)
    :m_value(value.m_value)
  {}
#endif
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline
  simd& operator=(storage_type const& value) {
    copy_from(value.data(), element_aligned_tag());
    return *this;
  }
  template <class Flags>
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd(T const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd(T const* ptr, int /*stride*/)
    : m_value(ptr[0])
  {}
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd operator*(simd const& other) const {
    return simd(m_value * other.m_value);
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd operator/(simd const& other) const {
    return simd(m_value / other.m_value);
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd operator+(simd const& other) const {
    return simd(m_value + other.m_value);
  }
#ifdef STK_VOLATILE_SIMD
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline void plus_equals(simd const volatile& other) volatile {
    m_value = m_value + other.m_value;
  }
#endif
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd operator-(simd const& other) const {
    return simd(m_value - other.m_value);
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd operator-() const {
    return simd(-m_value);
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE void copy_from(T const* ptr, element_aligned_tag) {
    m_value = *ptr;
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE void copy_to(T* ptr, element_aligned_tag) const {
    *ptr = m_value;
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE constexpr T get() const { return m_value; }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd_mask<T, simd_abi::scalar> operator<(simd const& other) const {
    return simd_mask<T, simd_abi::scalar>(m_value < other.m_value);
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd_mask<T, simd_abi::scalar> operator==(simd const& other) const {
    return simd_mask<T, simd_abi::scalar>(m_value == other.m_value);
  }
};

template <class T>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::scalar> abs(simd<T, simd_abi::scalar> const& a) {
  return simd<T, simd_abi::scalar>(std::abs(a.get()));
}

template <class T>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::scalar> sqrt(simd<T, simd_abi::scalar> const& a) {
  return simd<T, simd_abi::scalar>(std::sqrt(a.get()));
}

template <class T>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::scalar> cbrt(simd<T, simd_abi::scalar> const& a) {
  return simd<T, simd_abi::scalar>(std::cbrt(a.get()));
}

template <class T>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::scalar> exp(simd<T, simd_abi::scalar> const& a) {
  return simd<T, simd_abi::scalar>(std::exp(a.get()));
}

template <class T>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::scalar> fma(
    simd<T, simd_abi::scalar> const& a,
    simd<T, simd_abi::scalar> const& b,
    simd<T, simd_abi::scalar> const& c) {
  return simd<T, simd_abi::scalar>((a.get() * b.get()) + c.get());
}

template <class T>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::scalar> max(
    simd<T, simd_abi::scalar> const& a, simd<T, simd_abi::scalar> const& b) {
  return simd<T, simd_abi::scalar>(choose((a.get() < b.get()), b.get(), a.get()));
}

template <class T>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::scalar> min(
    simd<T, simd_abi::scalar> const& a, simd<T, simd_abi::scalar> const& b) {
  return simd<T, simd_abi::scalar>(choose((b.get() < a.get()), b.get(), a.get()));
}

template <class T>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, simd_abi::scalar> choose(
    simd_mask<T, simd_abi::scalar> const& a, simd<T, simd_abi::scalar> const& b, simd<T, simd_abi::scalar> const& c) {
  return simd<T, simd_abi::scalar>(choose(a.get(), b.get(), c.get()));
}

}

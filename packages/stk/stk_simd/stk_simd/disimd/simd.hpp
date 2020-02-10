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

#pragma once

#ifdef __SSE__
#include <xmmintrin.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#ifdef __AVX__
#include <immintrin.h>
#endif

#ifdef __AVX512F__
#include <immintrin.h>
#endif

#ifdef __ARM_NEON
#include <arm_neon.h>
#endif

#ifdef __VSX__
#include <altivec.h>
// undefine the really dangerous macros from this file
#undef vector
#undef pixel
#undef bool
#endif

#include <cmath>
#include <cstdint>

#ifdef __CUDACC__
#define SIMD_CUDA_ALWAYS_INLINE __forceinline__
#endif

#ifndef SIMD_ALWAYS_INLINE
#define SIMD_ALWAYS_INLINE [[gnu::always_inline]]
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

#ifndef SIMD_PRAGMA
#if defined(_OPENMP)
#define SIMD_PRAGMA _Pragma("omp simd")
#elif defined(__clang__)
#define SIMD_PRAGMA _Pragma("clang loop vectorize(enable)")
#elif defined(__GNUC__)
#define SIMD_PRAGMA _Pragma("GCC ivdep")
#else
#define SIMD_PRAGMA
#endif
#endif

#ifndef SIMD_NAMESPACE
#define SIMD_NAMESPACE simd
#endif

namespace SIMD_NAMESPACE {

template <class T, class Abi>
class simd;

template <class T, class Abi>
class simd_mask;

class element_aligned_tag {};

#ifndef SIMD_SCALAR_CHOOSE_DEFINED
template <class T>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE constexpr T const&
choose(bool a, T const& b, T const& c) {
  return a ? b : c;
}
#endif

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi>& operator+=(simd<T, Abi>& a, simd<T, Abi> const& b) {
  a = a + b;
  return a;
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi>& operator-=(simd<T, Abi>& a, simd<T, Abi> const& b) {
  a = a - b;
  return a;
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi>& operator*=(simd<T, Abi>& a, simd<T, Abi> const& b) {
  a = a * b;
  return a;
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi>& operator/=(simd<T, Abi>& a, simd<T, Abi> const& b) {
  a = a / b;
  return a;
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi> operator+(T const& a, simd<T, Abi> const& b) {
  return simd<T, Abi>(a) + b;
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi> operator+(simd<T, Abi> const& a, T const& b) {
  return a + simd<T, Abi>(b);
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi> operator-(T const& a, simd<T, Abi> const& b) {
  return simd<T, Abi>(a) - b;
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi> operator-(simd<T, Abi> const& a, T const& b) {
  return a - simd<T, Abi>(b);
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi> operator*(T const& a, simd<T, Abi> const& b) {
  return simd<T, Abi>(a) * b;
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi> operator*(simd<T, Abi> const& a, T const& b) {
  return a * simd<T, Abi>(b);
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi> operator/(T const& a, simd<T, Abi> const& b) {
  return simd<T, Abi>(a) / b;
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi> operator/(simd<T, Abi> const& a, T const& b) {
  return a / simd<T, Abi>(b);
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi> cbrt(simd<T, Abi> a) {
  T tmp[simd<T, Abi>::size()];
  a.copy_to(tmp, element_aligned_tag());
  for (int i = 0; i < simd<T, Abi>::size(); ++i) tmp[i] = std::cbrt(tmp[i]);
  a.copy_from(tmp, element_aligned_tag());
  return a;
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi> exp(simd<T, Abi> a) {
  T tmp[simd<T, Abi>::size()];
  a.copy_to(tmp, element_aligned_tag());
  for (int i = 0; i < simd<T, Abi>::size(); ++i) tmp[i] = std::exp(tmp[i]);
  a.copy_from(tmp, element_aligned_tag());
  return a;
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi> fma(simd<T, Abi> a, simd<T, Abi> const& b, simd<T, Abi> const& c) {
  T stack_a[simd<T, Abi>::size()];
  T stack_b[simd<T, Abi>::size()];
  a.copy_to(stack_a, element_aligned_tag());
  b.copy_to(stack_b, element_aligned_tag());
  for (int i = 0; i < simd<T, Abi>::size(); ++i) stack_a[i] *= stack_b[i];
  c.copy_to(stack_b, element_aligned_tag());
  for (int i = 0; i < simd<T, Abi>::size(); ++i) stack_a[i] += stack_b[i];
  a.copy_from(stack_a, element_aligned_tag());
  return a;
}

SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline bool all_of(bool a) { return a; }
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline bool any_of(bool a) { return a; }

template <class T, class Abi>
class simd_storage {
  T m_value[simd<T, Abi>::size()];
 public:
  using value_type = T;
  using simd_type = simd<T, Abi>;
  SIMD_ALWAYS_INLINE inline simd_storage() = default;
  SIMD_ALWAYS_INLINE inline static constexpr
  int size() { return simd<T, Abi>::size(); }
  SIMD_ALWAYS_INLINE inline
  simd_storage(simd<T, Abi> const& value) {
    value.copy_to(m_value, element_aligned_tag());
  }
  SIMD_ALWAYS_INLINE explicit inline
  simd_storage(T value)
    :simd_storage(simd<T, Abi>(value))
  {}
  SIMD_ALWAYS_INLINE inline
  simd_storage& operator=(simd<T, Abi> const& value) {
    value.copy_to(m_value, element_aligned_tag());
    return *this;
  }
  SIMD_ALWAYS_INLINE inline
  T const* data() const { return m_value; }
  SIMD_ALWAYS_INLINE inline
  T* data() { return m_value; }
  SIMD_ALWAYS_INLINE inline
  T const& operator[](int i) const { return m_value[i]; }
  SIMD_ALWAYS_INLINE inline
  T& operator[](int i) { return m_value[i]; }
};

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
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE static constexpr int size() { return 1; }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd(T value)
    :m_value(value)
  {}
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
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd(T const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd operator*(simd const& other) const {
    return simd(m_value * other.m_value);
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd operator/(simd const& other) const {
    return simd(m_value / other.m_value);
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd operator+(simd const& other) const {
    return simd(m_value + other.m_value);
  }
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
  using std::sqrt;
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

#if defined(__clang__)

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

#endif

/* Intel SVML disclaimer: cbrt, exp, etc. are not intrinsics, they are Intel-proprietary library functions
  https://stackoverflow.com/questions/36636159/where-is-clangs-mm256-pow-ps-intrinsic
  This is why the specializations that call these functions are protected with __INTEL_COMPILER.
 */

/* Intel FMA disclaimer: it is hard to detect FMA across compilers
   https://stackoverflow.com/questions/16348909/how-do-i-know-if-i-can-compile-with-fma-instruction-sets
   it seems like the best we can do is __FMA__ or __AVX2__, since MSVC doesn't define __FMA__
 */

#ifdef __SSE__

namespace simd_abi {

class sse {};

}

template <>
class simd_mask<float, simd_abi::sse> {
  __m128 m_value;
 public:
  using value_type = bool;
  using simd_type = simd<float, simd_abi::sse>;
  using abi_type = simd_abi::sse;
  SIMD_ALWAYS_INLINE inline simd_mask() = default;
  SIMD_ALWAYS_INLINE inline simd_mask(bool value)
    :m_value(_mm_castsi128_ps(_mm_set1_epi32(-int(value))))
  {}
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 4; }
  SIMD_ALWAYS_INLINE inline constexpr simd_mask(__m128 const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE constexpr __m128 get() const { return m_value; }
  SIMD_ALWAYS_INLINE simd_mask operator||(simd_mask const& other) const {
    return simd_mask(_mm_or_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE simd_mask operator&&(simd_mask const& other) const {
    return simd_mask(_mm_and_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE simd_mask operator!() const {
    return simd_mask(_mm_andnot_ps(m_value, simd_mask(true).get()));
  }
};

SIMD_ALWAYS_INLINE inline bool all_of(simd_mask<float, simd_abi::sse> const& a) {
  return _mm_movemask_ps(a.get()) == 0xF;
}

SIMD_ALWAYS_INLINE inline bool any_of(simd_mask<float, simd_abi::sse> const& a) {
  return _mm_movemask_ps(a.get()) != 0x0;
}

template <>
class simd<float, simd_abi::sse> {
  __m128 m_value;
 public:
  using value_type = float;
  using abi_type = simd_abi::sse;
  using mask_type = simd_mask<float, abi_type>;
  using storage_type = simd_storage<float, abi_type>;
  SIMD_ALWAYS_INLINE inline simd() = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 4; }
  SIMD_ALWAYS_INLINE inline simd(float value)
    :m_value(_mm_set1_ps(value))
  {}
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
  SIMD_ALWAYS_INLINE inline simd(float const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_ALWAYS_INLINE inline constexpr simd(__m128 const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline simd operator*(simd const& other) const {
    return simd(_mm_mul_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator/(simd const& other) const {
    return simd(_mm_div_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator+(simd const& other) const {
    return simd(_mm_add_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator-(simd const& other) const {
    return simd(_mm_sub_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator-() const {
    return simd(_mm_sub_ps(_mm_set1_ps(0.0), m_value));
  }
  SIMD_ALWAYS_INLINE void copy_from(float const* ptr, element_aligned_tag) {
    m_value = _mm_loadu_ps(ptr);
  }
  SIMD_ALWAYS_INLINE void copy_to(float* ptr, element_aligned_tag) const {
    _mm_storeu_ps(ptr, m_value);
  }
  SIMD_ALWAYS_INLINE constexpr __m128 get() const { return m_value; }
  SIMD_ALWAYS_INLINE simd_mask<float, simd_abi::sse> operator<(simd const& other) const {
    return simd_mask<float, simd_abi::sse>(_mm_cmplt_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE simd_mask<float, simd_abi::sse> operator==(simd const& other) const {
    return simd_mask<float, simd_abi::sse>(_mm_cmpeq_ps(m_value, other.m_value));
  }
};

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::sse> abs(simd<float, simd_abi::sse> const& a) {
  __m128 const sign_mask = _mm_set1_ps(-0.f);  // -0.f = 1 << 31
  return simd<float, simd_abi::sse>(_mm_andnot_ps(sign_mask, a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::sse> sqrt(simd<float, simd_abi::sse> const& a) {
  return simd<float, simd_abi::sse>(_mm_sqrt_ps(a.get()));
}

#ifdef __INTEL_COMPILER
SIMD_ALWAYS_INLINE inline simd<float, simd_abi::sse> cbrt(simd<float, simd_abi::sse> const& a) {
  return simd<float, simd_abi::sse>(_mm_cbrt_ps(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::sse> exp(simd<float, simd_abi::sse> const& a) {
  return simd<float, simd_abi::sse>(_mm_exp_ps(a.get()));
}
#endif

#if defined(__FMA__) || defined(__AVX2__)
SIMD_ALWAYS_INLINE inline simd<float, simd_abi::sse> fma(
    simd<float, simd_abi::sse> const& a,
    simd<float, simd_abi::sse> const& b,
    simd<float, simd_abi::sse> const& c) {
  return simd<float, simd_abi::sse>(_mm_fmadd_ps(a.get(), b.get(), c.get()));
}
#endif

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::sse> max(
    simd<float, simd_abi::sse> const& a, simd<float, simd_abi::sse> const& b) {
  return simd<float, simd_abi::sse>(_mm_max_ps(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::sse> min(
    simd<float, simd_abi::sse> const& a, simd<float, simd_abi::sse> const& b) {
  return simd<float, simd_abi::sse>(_mm_min_ps(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::sse> choose(
    simd_mask<float, simd_abi::sse> const& a, simd<float, simd_abi::sse> const& b, simd<float, simd_abi::sse> const& c) {
  return simd<float, simd_abi::sse>(_mm_add_ps(_mm_and_ps(a.get(), b.get()), _mm_andnot_ps(a.get(), c.get())));
}

#endif

#ifdef __SSE2__

template <>
class simd_mask<double, simd_abi::sse> {
  __m128d m_value;
 public:
  using value_type = bool;
  using simd_type = simd<double, simd_abi::sse>;
  using abi_type = simd_abi::sse;
  SIMD_ALWAYS_INLINE inline simd_mask() = default;
  SIMD_ALWAYS_INLINE inline simd_mask(bool value)
    :m_value(_mm_castsi128_pd(_mm_set1_epi64x(-std::int64_t(value))))
  {}
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 4; }
  SIMD_ALWAYS_INLINE inline constexpr simd_mask(__m128d const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline constexpr __m128d get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask operator||(simd_mask const& other) const {
    return simd_mask(_mm_or_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator&&(simd_mask const& other) const {
    return simd_mask(_mm_and_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator!() const {
    return simd_mask(_mm_andnot_pd(m_value, simd_mask(true).get()));
  }
};

SIMD_ALWAYS_INLINE inline bool all_of(simd_mask<double, simd_abi::sse> const& a) {
  return _mm_movemask_pd(a.get()) == 0x3;
}

SIMD_ALWAYS_INLINE inline bool any_of(simd_mask<double, simd_abi::sse> const& a) {
  return _mm_movemask_pd(a.get()) != 0x0;
}

template <>
class simd<double, simd_abi::sse> {
  __m128d m_value;
 public:
  using value_type = double;
  using abi_type = simd_abi::sse;
  using mask_type = simd_mask<double, abi_type>;
  using storage_type = simd_storage<double, abi_type>;
  SIMD_ALWAYS_INLINE inline simd() = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 2; }
  SIMD_ALWAYS_INLINE inline simd(double value)
    :m_value(_mm_set1_pd(value))
  {}
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
  SIMD_ALWAYS_INLINE inline simd(double const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_ALWAYS_INLINE inline constexpr simd(__m128d const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline simd operator*(simd const& other) const {
    return simd(_mm_mul_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator/(simd const& other) const {
    return simd(_mm_div_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator+(simd const& other) const {
    return simd(_mm_add_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator-(simd const& other) const {
    return simd(_mm_sub_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator-() const {
    return simd(_mm_sub_pd(_mm_set1_pd(0.0), m_value));
  }
  SIMD_ALWAYS_INLINE inline void copy_from(double const* ptr, element_aligned_tag) {
    m_value = _mm_loadu_pd(ptr);
  }
  SIMD_ALWAYS_INLINE inline void copy_to(double* ptr, element_aligned_tag) const {
    _mm_storeu_pd(ptr, m_value);
  }
  SIMD_ALWAYS_INLINE inline constexpr __m128d get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask<double, simd_abi::sse> operator<(simd const& other) const {
    return simd_mask<double, simd_abi::sse>(_mm_cmplt_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask<double, simd_abi::sse> operator==(simd const& other) const {
    return simd_mask<double, simd_abi::sse>(_mm_cmpeq_pd(m_value, other.m_value));
  }
};

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::sse> abs(simd<double, simd_abi::sse> const& a) {
  __m128d const sign_mask = _mm_set1_pd(-0.);  // -0. = 1 << 63
  return simd<double, simd_abi::sse>(_mm_andnot_pd(sign_mask, a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::sse> sqrt(simd<double, simd_abi::sse> const& a) {
  return simd<double, simd_abi::sse>(_mm_sqrt_pd(a.get()));
}

#ifdef __INTEL_COMPILER
SIMD_ALWAYS_INLINE inline simd<double, simd_abi::sse> cbrt(simd<double, simd_abi::sse> const& a) {
  return simd<double, simd_abi::sse>(_mm_cbrt_pd(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::sse> exp(simd<double, simd_abi::sse> const& a) {
  return simd<double, simd_abi::sse>(_mm_exp_pd(a.get()));
}
#endif

#if defined(__FMA__) || defined(__AVX2__)
SIMD_ALWAYS_INLINE inline simd<double, simd_abi::sse> fma(
    simd<double, simd_abi::sse> const& a,
    simd<double, simd_abi::sse> const& b,
    simd<double, simd_abi::sse> const& c) {
  return simd<double, simd_abi::sse>(_mm_fmadd_pd(a.get(), b.get(), c.get()));
}
#endif

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::sse> max(
    simd<double, simd_abi::sse> const& a, simd<double, simd_abi::sse> const& b) {
  return simd<double, simd_abi::sse>(_mm_max_pd(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::sse> min(
    simd<double, simd_abi::sse> const& a, simd<double, simd_abi::sse> const& b) {
  return simd<double, simd_abi::sse>(_mm_min_pd(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::sse> choose(
    simd_mask<double, simd_abi::sse> const& a, simd<double, simd_abi::sse> const& b, simd<double, simd_abi::sse> const& c) {
  return simd<double, simd_abi::sse>(
      _mm_add_pd(
        _mm_and_pd(a.get(), b.get()),
        _mm_andnot_pd(a.get(), c.get())));
}

#endif

#ifdef __AVX__

namespace simd_abi {

class avx {};

}

template <>
class simd_mask<float, simd_abi::avx> {
  __m256 m_value;
 public:
  using value_type = bool;
  using simd_type = simd<float, simd_abi::avx>;
  using abi_type = simd_abi::avx;
  SIMD_ALWAYS_INLINE inline simd_mask() = default;
  SIMD_ALWAYS_INLINE inline simd_mask(bool value) {
    m_value = _mm256_castsi256_ps(_mm256_set1_epi32(-int(value)));
  }
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 8; }
  SIMD_ALWAYS_INLINE inline constexpr simd_mask(__m256 const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline constexpr __m256 get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask operator||(simd_mask const& other) const {
    return simd_mask(_mm256_or_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator&&(simd_mask const& other) const {
    return simd_mask(_mm256_and_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator!() const {
    return simd_mask(_mm256_andnot_ps(m_value, simd_mask(true).get()));
  }
};

SIMD_ALWAYS_INLINE inline bool all_of(simd_mask<float, simd_abi::avx> const& a) {
  return _mm256_testc_ps(a.get(), simd_mask<float, simd_abi::avx>(true).get());
}

SIMD_ALWAYS_INLINE inline bool any_of(simd_mask<float, simd_abi::avx> const& a) {
  return !_mm256_testc_ps(simd_mask<float, simd_abi::avx>(false).get(), a.get());
}

template <>
class simd<float, simd_abi::avx> {
  __m256 m_value;
 public:
  using value_type = float;
  using abi_type = simd_abi::avx;
  using mask_type = simd_mask<float, abi_type>;
  using storage_type = simd_storage<float, abi_type>;
  SIMD_ALWAYS_INLINE inline simd() = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 8; }
  SIMD_ALWAYS_INLINE inline simd(float value)
    :m_value(_mm256_set1_ps(value))
  {}
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
  SIMD_ALWAYS_INLINE inline simd(float const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_ALWAYS_INLINE inline constexpr simd(__m256 const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline simd operator*(simd const& other) const {
    return simd(_mm256_mul_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator/(simd const& other) const {
    return simd(_mm256_div_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator+(simd const& other) const {
    return simd(_mm256_add_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator-(simd const& other) const {
    return simd(_mm256_sub_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd operator-() const {
    return simd(_mm256_sub_ps(_mm256_set1_ps(0.0), m_value));
  }
  SIMD_ALWAYS_INLINE inline void copy_from(float const* ptr, element_aligned_tag) {
    m_value = _mm256_loadu_ps(ptr);
  }
  SIMD_ALWAYS_INLINE inline void copy_to(float* ptr, element_aligned_tag) const {
    _mm256_storeu_ps(ptr, m_value);
  }
  SIMD_ALWAYS_INLINE inline constexpr __m256 get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask<float, simd_abi::avx> operator<(simd const& other) const {
    return simd_mask<float, simd_abi::avx>(_mm256_cmp_ps(m_value, other.m_value, _CMP_LT_OS));
  }
  SIMD_ALWAYS_INLINE inline simd_mask<float, simd_abi::avx> operator==(simd const& other) const {
    return simd_mask<float, simd_abi::avx>(_mm256_cmp_ps(m_value, other.m_value, _CMP_EQ_OS));
  }
};

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx> abs(simd<float, simd_abi::avx> const& a) {
  __m256 sign_mask = _mm256_set1_ps(-0.f);  // -0.f = 1 << 31
  return simd<float, simd_abi::avx>(_mm256_andnot_ps(sign_mask, a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx> sqrt(simd<float, simd_abi::avx> const& a) {
  return simd<float, simd_abi::avx>(_mm256_sqrt_ps(a.get()));
}

#ifdef __INTEL_COMPILER
SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx> cbrt(simd<float, simd_abi::avx> const& a) {
  return simd<float, simd_abi::avx>(_mm256_cbrt_ps(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx> exp(simd<float, simd_abi::avx> const& a) {
  return simd<float, simd_abi::avx>(_mm256_exp_ps(a.get()));
}
#endif

#if defined(__FMA__) || defined(__AVX2__)
SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx> fma(
    simd<float, simd_abi::avx> const& a,
    simd<float, simd_abi::avx> const& b,
    simd<float, simd_abi::avx> const& c) {
  return simd<float, simd_abi::avx>(_mm256_fmadd_ps(a.get(), b.get(), c.get()));
}
#endif

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx> max(
    simd<float, simd_abi::avx> const& a, simd<float, simd_abi::avx> const& b) {
  return simd<float, simd_abi::avx>(_mm256_max_ps(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx> min(
    simd<float, simd_abi::avx> const& a, simd<float, simd_abi::avx> const& b) {
  return simd<float, simd_abi::avx>(_mm256_min_ps(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx> choose(
    simd_mask<float, simd_abi::avx> const& a, simd<float, simd_abi::avx> const& b, simd<float, simd_abi::avx> const& c) {
  return simd<float, simd_abi::avx>(_mm256_blendv_ps(c.get(), b.get(), a.get()));
}

template <>
class simd_mask<double, simd_abi::avx> {
  __m256d m_value;
 public:
  using value_type = bool;
  using simd_type = simd<double, simd_abi::avx>;
  using abi_type = simd_abi::avx;
  SIMD_ALWAYS_INLINE inline simd_mask() = default;
  SIMD_ALWAYS_INLINE inline simd_mask(bool value) {
    m_value = _mm256_castsi256_pd(_mm256_set1_epi64x(-std::int64_t(value)));
  }
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 4; }
  SIMD_ALWAYS_INLINE inline constexpr simd_mask(__m256d const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline constexpr __m256d get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask operator||(simd_mask const& other) const {
    return simd_mask(_mm256_or_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator&&(simd_mask const& other) const {
    return simd_mask(_mm256_and_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator!() const {
    return simd_mask(_mm256_andnot_pd(m_value, simd_mask(true).get()));
  }
};

SIMD_ALWAYS_INLINE inline bool all_of(simd_mask<double, simd_abi::avx> const& a) {
  return _mm256_testc_pd(a.get(),
      simd_mask<double, simd_abi::avx>(true).get());
}

SIMD_ALWAYS_INLINE inline bool any_of(simd_mask<double, simd_abi::avx> const& a) {
  return !_mm256_testc_pd(
      simd_mask<double, simd_abi::avx>(false).get(), a.get());
}

template <>
class simd<double, simd_abi::avx> {
  __m256d m_value;
 public:
  using value_type = double;
  using abi_type = simd_abi::avx;
  using mask_type = simd_mask<double, abi_type>;
  using storage_type = simd_storage<double, abi_type>;
  SIMD_ALWAYS_INLINE inline simd() = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 4; }
  SIMD_ALWAYS_INLINE inline simd(double value)
    :m_value(_mm256_set1_pd(value))
  {}
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
  SIMD_ALWAYS_INLINE inline simd(double const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_ALWAYS_INLINE inline constexpr simd(__m256d const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline simd operator*(simd const& other) const {
    return simd(_mm256_mul_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator/(simd const& other) const {
    return simd(_mm256_div_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator+(simd const& other) const {
    return simd(_mm256_add_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator-(simd const& other) const {
    return simd(_mm256_sub_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd operator-() const {
    return simd(_mm256_sub_pd(_mm256_set1_pd(0.0), m_value));
  }
  SIMD_ALWAYS_INLINE inline void copy_from(double const* ptr, element_aligned_tag) {
    m_value = _mm256_loadu_pd(ptr);
  }
  SIMD_ALWAYS_INLINE inline void copy_to(double* ptr, element_aligned_tag) const {
    _mm256_storeu_pd(ptr, m_value);
  }
  SIMD_ALWAYS_INLINE inline constexpr __m256d get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask<double, simd_abi::avx> operator<(simd const& other) const {
    return simd_mask<double, simd_abi::avx>(_mm256_cmp_pd(m_value, other.m_value, _CMP_LT_OS));
  }
  SIMD_ALWAYS_INLINE inline simd_mask<double, simd_abi::avx> operator==(simd const& other) const {
    return simd_mask<double, simd_abi::avx>(_mm256_cmp_pd(m_value, other.m_value, _CMP_EQ_OS));
  }
};

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx> abs(simd<double, simd_abi::avx> const& a) {
  __m256d const sign_mask = _mm256_set1_pd(-0.f);  // -0.f = 1 << 31
  return simd<double, simd_abi::avx>(_mm256_andnot_pd(sign_mask, a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx> sqrt(simd<double, simd_abi::avx> const& a) {
  return simd<double, simd_abi::avx>(_mm256_sqrt_pd(a.get()));
}

#ifdef __INTEL_COMPILER
SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx> cbrt(simd<double, simd_abi::avx> const& a) {
  return simd<double, simd_abi::avx>(_mm256_cbrt_pd(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx> exp(simd<double, simd_abi::avx> const& a) {
  return simd<double, simd_abi::avx>(_mm256_exp_pd(a.get()));
}
#endif

#if defined(__FMA__) || defined(__AVX2__)
SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx> fma(
    simd<double, simd_abi::avx> const& a,
    simd<double, simd_abi::avx> const& b,
    simd<double, simd_abi::avx> const& c) {
  return simd<double, simd_abi::avx>(_mm256_fmadd_pd(a.get(), b.get(), c.get()));
}
#endif

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx> max(
    simd<double, simd_abi::avx> const& a, simd<double, simd_abi::avx> const& b) {
  return simd<double, simd_abi::avx>(_mm256_max_pd(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx> min(
    simd<double, simd_abi::avx> const& a, simd<double, simd_abi::avx> const& b) {
  return simd<double, simd_abi::avx>(_mm256_min_pd(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx> choose(
    simd_mask<double, simd_abi::avx> const& a, simd<double, simd_abi::avx> const& b, simd<double, simd_abi::avx> const& c) {
  return simd<double, simd_abi::avx>(_mm256_blendv_pd(c.get(), b.get(), a.get()));
}

#endif

#ifdef __AVX512F__

namespace simd_abi {

class avx512 {};

}

template <>
class simd_mask<float, simd_abi::avx512> {
  __mmask16 m_value;
 public:
  using value_type = bool;
  using simd_type = simd<float, simd_abi::avx512>;
  using abi_type = simd_abi::avx512;
  SIMD_ALWAYS_INLINE inline simd_mask() = default;
  SIMD_ALWAYS_INLINE inline simd_mask(bool value)
    :m_value(-std::int16_t(value))
  {}
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 16; }
  SIMD_ALWAYS_INLINE inline constexpr simd_mask(__mmask16 const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline constexpr __mmask16 get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask operator||(simd_mask const& other) const {
    return simd_mask(_kor_mask16(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator&&(simd_mask const& other) const {
    return simd_mask(_kand_mask16(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator!() const {
    return simd_mask(_knot_mask16(m_value));
  }
};

SIMD_ALWAYS_INLINE inline bool all_of(simd_mask<float, simd_abi::avx512> const& a) {
  return _ktestc_mask16_u8(a.get(),
      simd_mask<float, simd_abi::avx512>(true).get());
}

SIMD_ALWAYS_INLINE inline bool any_of(simd_mask<float, simd_abi::avx512> const& a) {
  return !_ktestc_mask16_u8(
      simd_mask<float, simd_abi::avx512>(false).get(), a.get());
}

template <>
class simd<float, simd_abi::avx512> {
  __m512 m_value;
 public:
  SIMD_ALWAYS_INLINE simd() = default;
  using value_type = float;
  using abi_type = simd_abi::avx512;
  using mask_type = simd_mask<float, abi_type>;
  using storage_type = simd_storage<float, abi_type>;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 16; }
  SIMD_ALWAYS_INLINE inline simd(float value)
    :m_value(_mm512_set1_ps(value))
  {}
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
  SIMD_ALWAYS_INLINE inline simd(float const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_ALWAYS_INLINE inline constexpr simd(__m512 const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline simd operator*(simd const& other) const {
    return simd(_mm512_mul_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator/(simd const& other) const {
    return simd(_mm512_div_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator+(simd const& other) const {
    return simd(_mm512_add_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator-(simd const& other) const {
    return simd(_mm512_sub_ps(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd operator-() const {
    return simd(_mm512_sub_ps(_mm512_set1_ps(0.0), m_value));
  }
  SIMD_ALWAYS_INLINE inline void copy_from(float const* ptr, element_aligned_tag) {
    m_value = _mm512_loadu_ps(ptr);
  }
  SIMD_ALWAYS_INLINE inline void copy_to(float* ptr, element_aligned_tag) const {
    _mm512_storeu_ps(ptr, m_value);
  }
  SIMD_ALWAYS_INLINE inline constexpr __m512 get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask<float, simd_abi::avx512> operator<(simd const& other) const {
    return simd_mask<float, simd_abi::avx512>(_mm512_cmp_ps_mask(m_value, other.m_value, _CMP_LT_OS));
  }
  SIMD_ALWAYS_INLINE inline simd_mask<float, simd_abi::avx512> operator==(simd const& other) const {
    return simd_mask<float, simd_abi::avx512>(_mm512_cmp_ps_mask(m_value, other.m_value, _CMP_EQ_OS));
  }
};

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx512> abs(simd<float, simd_abi::avx512> const& a) {
  __m512 const rhs = a.get();
  return reinterpret_cast<__m512>(_mm512_and_epi32(reinterpret_cast<__m512i>(rhs), _mm512_set1_epi32(0x7fffffff)));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx512> sqrt(simd<float, simd_abi::avx512> const& a) {
  return simd<float, simd_abi::avx512>(_mm512_sqrt_ps(a.get()));
}

#ifdef __INTEL_COMPILER
SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx512> cbrt(simd<float, simd_abi::avx512> const& a) {
  return simd<float, simd_abi::avx512>(_mm512_cbrt_ps(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx512> exp(simd<float, simd_abi::avx512> const& a) {
  return simd<float, simd_abi::avx512>(_mm512_exp_ps(a.get()));
}
#endif

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx512> fma(
    simd<float, simd_abi::avx512> const& a,
    simd<float, simd_abi::avx512> const& b,
    simd<float, simd_abi::avx512> const& c) {
  return simd<float, simd_abi::avx512>(_mm512_fmadd_ps(a.get(), b.get(), c.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx512> max(
    simd<float, simd_abi::avx512> const& a, simd<float, simd_abi::avx512> const& b) {
  return simd<float, simd_abi::avx512>(_mm512_max_ps(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx512> min(
    simd<float, simd_abi::avx512> const& a, simd<float, simd_abi::avx512> const& b) {
  return simd<float, simd_abi::avx512>(_mm512_min_ps(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx512> choose(
    simd_mask<float, simd_abi::avx512> const& a, simd<float, simd_abi::avx512> const& b, simd<float, simd_abi::avx512> const& c) {
  return simd<float, simd_abi::avx512>(_mm512_mask_blend_ps(a.get(), c.get(), b.get()));
}

template <>
class simd_mask<double, simd_abi::avx512> {
  __mmask8 m_value;
 public:
  using value_type = bool;
  SIMD_ALWAYS_INLINE inline simd_mask() = default;
  SIMD_ALWAYS_INLINE inline simd_mask(bool value)
    :m_value(-std::int16_t(value))
  {}
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 8; }
  SIMD_ALWAYS_INLINE inline constexpr simd_mask(__mmask8 const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline constexpr __mmask8 get() const { return m_value; }
  SIMD_ALWAYS_INLINE simd_mask operator||(simd_mask const& other) const {
    return simd_mask(_kor_mask8(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE simd_mask operator&&(simd_mask const& other) const {
    return simd_mask(_kand_mask8(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE simd_mask operator!() const {
    return simd_mask(_knot_mask8(m_value));
  }
};

SIMD_ALWAYS_INLINE inline bool all_of(simd_mask<double, simd_abi::avx512> const& a) {
  return _ktestc_mask8_u8(a.get(),
      simd_mask<double, simd_abi::avx512>(true).get());
}

SIMD_ALWAYS_INLINE inline bool any_of(simd_mask<double, simd_abi::avx512> const& a) {
  return !_ktestc_mask8_u8(
      simd_mask<double, simd_abi::avx512>(false).get(), a.get());
}

template <>
class simd<double, simd_abi::avx512> {
  __m512d m_value;
 public:
  using value_type = double;
  using abi_type = simd_abi::avx512;
  using mask_type = simd_mask<double, abi_type>;
  using storage_type = simd_storage<double, abi_type>;
  SIMD_ALWAYS_INLINE inline simd() = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 8; }
  SIMD_ALWAYS_INLINE inline simd(double value)
    :m_value(_mm512_set1_pd(value))
  {}
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
  SIMD_ALWAYS_INLINE inline simd(double const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_ALWAYS_INLINE inline constexpr simd(__m512d const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline simd operator*(simd const& other) const {
    return simd(_mm512_mul_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator/(simd const& other) const {
    return simd(_mm512_div_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator+(simd const& other) const {
    return simd(_mm512_add_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator-(simd const& other) const {
    return simd(_mm512_sub_pd(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd operator-() const {
    return simd(_mm512_sub_pd(_mm512_set1_pd(0.0), m_value));
  }
  SIMD_ALWAYS_INLINE inline void copy_from(double const* ptr, element_aligned_tag) {
    m_value = _mm512_loadu_pd(ptr);
  }
  SIMD_ALWAYS_INLINE inline void copy_to(double* ptr, element_aligned_tag) const {
    _mm512_storeu_pd(ptr, m_value);
  }
  SIMD_ALWAYS_INLINE inline constexpr __m512d get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask<double, simd_abi::avx512> operator<(simd const& other) const {
    return simd_mask<double, simd_abi::avx512>(_mm512_cmp_pd_mask(m_value, other.m_value, _CMP_LT_OS));
  }
  SIMD_ALWAYS_INLINE inline simd_mask<double, simd_abi::avx512> operator==(simd const& other) const {
    return simd_mask<double, simd_abi::avx512>(_mm512_cmp_pd_mask(m_value, other.m_value, _CMP_EQ_OS));
  }
};

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx512> abs(simd<double, simd_abi::avx512> const& a) {
  __m512d const rhs = a.get();
  return reinterpret_cast<__m512d>(_mm512_and_epi64(_mm512_set1_epi64(0x7FFFFFFFFFFFFFFF),
        reinterpret_cast<__m512i>(rhs)));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx512> sqrt(simd<double, simd_abi::avx512> const& a) {
  return simd<double, simd_abi::avx512>(_mm512_sqrt_pd(a.get()));
}

#ifdef __INTEL_COMPILER
SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx512> cbrt(simd<double, simd_abi::avx512> const& a) {
  return simd<double, simd_abi::avx512>(_mm512_cbrt_pd(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx512> exp(simd<double, simd_abi::avx512> const& a) {
  return simd<double, simd_abi::avx512>(_mm512_exp_pd(a.get()));
}
#endif

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx512> fma(
    simd<double, simd_abi::avx512> const& a,
    simd<double, simd_abi::avx512> const& b,
    simd<double, simd_abi::avx512> const& c) {
  return simd<double, simd_abi::avx512>(_mm512_fmadd_pd(a.get(), b.get(), c.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx512> max(
    simd<double, simd_abi::avx512> const& a, simd<double, simd_abi::avx512> const& b) {
  return simd<double, simd_abi::avx512>(_mm512_max_pd(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx512> min(
    simd<double, simd_abi::avx512> const& a, simd<double, simd_abi::avx512> const& b) {
  return simd<double, simd_abi::avx512>(_mm512_min_pd(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx512> choose(
    simd_mask<double, simd_abi::avx512> const& a, simd<double, simd_abi::avx512> const& b, simd<double, simd_abi::avx512> const& c) {
  return simd<double, simd_abi::avx512>(_mm512_mask_blend_pd(a.get(), c.get(), b.get()));
}

#endif

#ifdef __ARM_NEON

namespace simd_abi {

class neon {};

}

template <>
class simd_mask<float, simd_abi::neon> {
  uint32x4_t m_value;
 public:
  using value_type = bool;
  using simd_type = simd_mask<float, simd_abi::neon>;
  using abi_type = simd_abi::neon;
  SIMD_ALWAYS_INLINE inline simd_mask() = default;
  SIMD_ALWAYS_INLINE inline simd_mask(bool value)
    :m_value(vreinterpretq_u32_s32(vdupq_n_s32(-int(value))))
  {}
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 4; }
  SIMD_ALWAYS_INLINE inline constexpr simd_mask(uint32x4_t const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline constexpr uint32x4_t get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask operator||(simd_mask const& other) const {
    return simd_mask(vorrq_u32(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator&&(simd_mask const& other) const {
    return simd_mask(vandq_u32(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator!() const {
    return simd_mask(vmvnq_u32(m_value));
  }
};

SIMD_ALWAYS_INLINE inline bool all_of(simd_mask<float, simd_abi::neon> const& a) {
  return vminvq_u32(a.get()) == std::uint32_t(-std::int32_t(1));
}

SIMD_ALWAYS_INLINE inline bool any_of(simd_mask<float, simd_abi::neon> const& a) {
  return vmaxvq_u32(a.get()) == std::uint32_t(-std::int32_t(1));
}

template <>
class simd<float, simd_abi::neon> {
  float32x4_t m_value;
 public:
  using value_type = float;
  using abi_type = simd_abi::neon;
  using mask_type = simd_mask<float, abi_type>;
  using storage_type = simd_storage<float, abi_type>;
  SIMD_ALWAYS_INLINE inline simd() = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 4; }
  SIMD_ALWAYS_INLINE inline simd(float value)
    :m_value(vdupq_n_f32(value))
  {}
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
  SIMD_ALWAYS_INLINE inline simd(float const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_ALWAYS_INLINE inline constexpr simd(float32x4_t const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline simd operator*(simd const& other) const {
    return simd(vmulq_f32(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator/(simd const& other) const {
    return simd(vdivq_f32(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator+(simd const& other) const {
    return simd(vaddq_f32(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator-(simd const& other) const {
    return simd(vsubq_f32(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator-() const {
    return simd(vnegq_f32(m_value));
  }
  SIMD_ALWAYS_INLINE inline void copy_from(float const* ptr, element_aligned_tag) {
    m_value = vld1q_f32(ptr);
  }
  SIMD_ALWAYS_INLINE inline void copy_to(float* ptr, element_aligned_tag) const {
    vst1q_f32(ptr, m_value);
  }
  SIMD_ALWAYS_INLINE inline constexpr float32x4_t get() const { return m_value; }
  SIMD_ALWAYS_INLINE simd_mask<float, simd_abi::neon> operator<(simd const& other) const {
    return simd_mask<float, simd_abi::neon>(vcltq_f32(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask<float, simd_abi::neon> operator==(simd const& other) const {
    return simd_mask<float, simd_abi::neon>(vceqq_f32(m_value, other.m_value));
  }
};

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::neon> abs(simd<float, simd_abi::neon> const& a) {
  return simd<float, simd_abi::neon>(vabsq_f32(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::neon> sqrt(simd<float, simd_abi::neon> const& a) {
  return simd<float, simd_abi::neon>(vsqrtq_f32(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::neon> fma(
    simd<float, simd_abi::neon> const& a,
    simd<float, simd_abi::neon> const& b,
    simd<float, simd_abi::neon> const& c) {
  return simd<float, simd_abi::neon>(vfmaq_f32(c.get(), b.get(), a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::neon> max(
    simd<float, simd_abi::neon> const& a, simd<float, simd_abi::neon> const& b) {
  return simd<float, simd_abi::neon>(vmaxq_f32(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::neon> min(
    simd<float, simd_abi::neon> const& a, simd<float, simd_abi::neon> const& b) {
  return simd<float, simd_abi::neon>(vminq_f32(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::neon> choose(
    simd_mask<float, simd_abi::neon> const& a, simd<float, simd_abi::neon> const& b, simd<float, simd_abi::neon> const& c) {
  return simd<float, simd_abi::neon>(
    vreinterpretq_f32_u32(
      vbslq_u32(
        a.get(),
        vreinterpretq_u32_f32(b.get()),
        vreinterpretq_u32_f32(c.get()))));
}

template <>
class simd_mask<double, simd_abi::neon> {
  uint64x2_t m_value;
 public:
  using value_type = bool;
  using simd_type = simd<double, simd_abi::neon>;
  using abi_type = simd_abi::neon;
  SIMD_ALWAYS_INLINE inline simd_mask() = default;
  SIMD_ALWAYS_INLINE inline simd_mask(bool value)
    :m_value(vreinterpretq_u64_s64(vdupq_n_s64(-std::int64_t(value))))
  {}
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 4; }
  SIMD_ALWAYS_INLINE inline constexpr simd_mask(uint64x2_t const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline constexpr uint64x2_t get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask operator||(simd_mask const& other) const {
    return simd_mask(vorrq_u64(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator&&(simd_mask const& other) const {
    return simd_mask(vandq_u64(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator!() const {
    return simd_mask(vreinterpretq_u64_u32(vmvnq_u32(vreinterpretq_u32_u64(m_value))));
  }
};

SIMD_ALWAYS_INLINE inline bool all_of(simd_mask<double, simd_abi::neon> const& a) {
  return all_of(simd_mask<float, simd_abi::neon>(vreinterpretq_u32_u64(a.get())));
}

SIMD_ALWAYS_INLINE inline bool any_of(simd_mask<double, simd_abi::neon> const& a) {
  return any_of(simd_mask<float, simd_abi::neon>(vreinterpretq_u32_u64(a.get())));
}

template <>
class simd<double, simd_abi::neon> {
  float64x2_t m_value;
 public:
  using value_type = double;
  using abi_type = simd_abi::neon;
  using mask_type = simd_mask<double, abi_type>;
  using storage_type = simd_storage<double, abi_type>;
  SIMD_ALWAYS_INLINE inline simd() = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 2; }
  SIMD_ALWAYS_INLINE inline simd(double value)
    :m_value(vdupq_n_f64(value))
  {}
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
  SIMD_ALWAYS_INLINE inline simd(double const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_ALWAYS_INLINE inline constexpr simd(float64x2_t const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline simd operator*(simd const& other) const {
    return simd(vmulq_f64(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator/(simd const& other) const {
    return simd(vdivq_f64(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator+(simd const& other) const {
    return simd(vaddq_f64(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator-(simd const& other) const {
    return simd(vsubq_f64(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator-() const {
    return simd(vnegq_f64(m_value));
  }
  SIMD_ALWAYS_INLINE inline void copy_from(double const* ptr, element_aligned_tag) {
    m_value = vld1q_f64(ptr);
  }
  SIMD_ALWAYS_INLINE inline void copy_to(double* ptr, element_aligned_tag) const {
    vst1q_f64(ptr, m_value);
  }
  SIMD_ALWAYS_INLINE inline constexpr float64x2_t get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask<double, simd_abi::neon> operator<(simd const& other) const {
    return simd_mask<double, simd_abi::neon>(vcltq_f64(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask<double, simd_abi::neon> operator==(simd const& other) const {
    return simd_mask<double, simd_abi::neon>(vceqq_f64(m_value, other.m_value));
  }
};

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::neon> abs(simd<double, simd_abi::neon> const& a) {
  return simd<double, simd_abi::neon>(vabsq_f64(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::neon> sqrt(simd<double, simd_abi::neon> const& a) {
  return simd<double, simd_abi::neon>(vsqrtq_f64(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::neon> fma(
    simd<double, simd_abi::neon> const& a,
    simd<double, simd_abi::neon> const& b,
    simd<double, simd_abi::neon> const& c) {
  return simd<double, simd_abi::neon>(vfmaq_f64(c.get(), b.get(), a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::neon> max(
    simd<double, simd_abi::neon> const& a, simd<double, simd_abi::neon> const& b) {
  return simd<double, simd_abi::neon>(vmaxq_f64(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::neon> min(
    simd<double, simd_abi::neon> const& a, simd<double, simd_abi::neon> const& b) {
  return simd<double, simd_abi::neon>(vminq_f64(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::neon> choose(
    simd_mask<double, simd_abi::neon> const& a, simd<double, simd_abi::neon> const& b, simd<double, simd_abi::neon> const& c) {
  return simd<double, simd_abi::neon>(
    vreinterpretq_f64_u64(
      vbslq_u64(
        a.get(),
        vreinterpretq_u64_f64(b.get()),
        vreinterpretq_u64_f64(c.get()))));
}

#endif

#if defined(__VSX__) && (!defined(__CUDACC__))

namespace simd_abi {

class vsx {};

}

template <>
class simd_mask<float, simd_abi::vsx> {
  __vector __bool int m_value;
 public:
  using value_type = bool;
  using simd_type = simd<float, simd_abi::vsx>;
  using abi_type = simd_abi::vsx;
  SIMD_ALWAYS_INLINE inline simd_mask() = default;
  SIMD_ALWAYS_INLINE inline simd_mask(bool value)
    :m_value{value, value, value, value}
  {}
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 4; }
  SIMD_ALWAYS_INLINE inline constexpr simd_mask(__vector __bool int const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline constexpr __vector __bool int get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask operator||(simd_mask const& other) const {
    return simd_mask(vec_or(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator&&(simd_mask const& other) const {
    return simd_mask(vec_and(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator!() const {
    return simd_mask(vec_nand(m_value, simd_mask(true).get()));
  }
};

SIMD_ALWAYS_INLINE inline bool all_of(simd_mask<float, simd_abi::vsx> const& a) {
  auto const true_value = simd_mask<float, simd_abi::vsx>(true).get();
  return vec_all_eq(a.get(), true_value);
}

SIMD_ALWAYS_INLINE inline bool any_of(simd_mask<float, simd_abi::vsx> const& a) {
  auto const true_value = simd_mask<float, simd_abi::vsx>(true).get();
  return vec_any_eq(a.get(), true_value);
}

template <>
class simd<float, simd_abi::vsx> {
  __vector float m_value;
 public:
  using value_type = float;
  using abi_type = simd_abi::vsx;
  using mask_type = simd_mask<float, abi_type>;
  using storage_type = simd_storage<float, abi_type>;
  SIMD_ALWAYS_INLINE inline simd() = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 4; }
  SIMD_ALWAYS_INLINE inline simd(float value)
    :m_value(vec_splats(value))
  {}
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
  SIMD_ALWAYS_INLINE inline simd(float const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_ALWAYS_INLINE inline constexpr simd(__vector float const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE simd operator*(simd const& other) const {
    return simd(vec_mul(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE simd operator/(simd const& other) const {
    return simd(vec_div(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE simd operator+(simd const& other) const {
    return simd(vec_add(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE simd operator-(simd const& other) const {
    return simd(vec_sub(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE simd operator-() const {
    // return simd(vec_neg(m_value)); some GCC versions dont have this
    return simd(0.0) - (*this);
  }
  SIMD_ALWAYS_INLINE void copy_from(float const* ptr, element_aligned_tag) {
    m_value = vec_vsx_ld(0, ptr);
  }
  SIMD_ALWAYS_INLINE void copy_to(float* ptr, element_aligned_tag) const {
    vec_vsx_st(m_value, 0, ptr);
  }
  SIMD_ALWAYS_INLINE constexpr __vector float get() const { return m_value; }
  SIMD_ALWAYS_INLINE simd_mask<float, simd_abi::vsx> operator<(simd const& other) const {
    return simd_mask<float, simd_abi::vsx>(vec_cmplt(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE simd_mask<float, simd_abi::vsx> operator==(simd const& other) const {
    return simd_mask<float, simd_abi::vsx>(vec_cmpeq(m_value, other.m_value));
  }
};

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::vsx> abs(simd<float, simd_abi::vsx> const& a) {
  return simd<float, simd_abi::vsx>(vec_abs(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::vsx> sqrt(simd<float, simd_abi::vsx> const& a) {
  return simd<float, simd_abi::vsx>(vec_sqrt(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::vsx> fma(
    simd<float, simd_abi::vsx> const& a,
    simd<float, simd_abi::vsx> const& b,
    simd<float, simd_abi::vsx> const& c) {
  return simd<float, simd_abi::vsx>(vec_madd(a.get(), b.get(), c.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::vsx> max(
    simd<float, simd_abi::vsx> const& a, simd<float, simd_abi::vsx> const& b) {
  return simd<float, simd_abi::vsx>(vec_max(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::vsx> min(
    simd<float, simd_abi::vsx> const& a, simd<float, simd_abi::vsx> const& b) {
  return simd<float, simd_abi::vsx>(vec_min(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::vsx> choose(
    simd_mask<float, simd_abi::vsx> const& a, simd<float, simd_abi::vsx> const& b, simd<float, simd_abi::vsx> const& c) {
  return simd<float, simd_abi::vsx>(vec_sel(c.get(), b.get(), a.get()));
}

template <>
class simd_mask<double, simd_abi::vsx> {
  /* Note: ideally, __vector __bool long would be the thing to use here.
   * however, GCC is missing key functions like vec_and for __vector __bool long.
   * that is why we use __vector unsigned long instead and convert back and forth
   */
  using ideal_type = __vector __bool long;
  using supported_type = __vector unsigned long;
  supported_type m_value;
  using ulong_t = unsigned long;
 public:
  using value_type = bool;
  using simd_type = simd_mask<double, simd_abi::vsx>;
  using abi_type = simd_abi::vsx;
  SIMD_ALWAYS_INLINE inline simd_mask() = default;
  SIMD_ALWAYS_INLINE inline simd_mask(bool value)
    :m_value{ulong_t(-long(value)), ulong_t(-long(value))}
  {}
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 2; }
  SIMD_ALWAYS_INLINE inline constexpr simd_mask(ideal_type const& value_in)
    :m_value(supported_type(value_in))
  {}
  SIMD_ALWAYS_INLINE inline constexpr simd_mask(supported_type const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline constexpr supported_type get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask operator||(simd_mask const& other) const {
    return simd_mask(vec_or(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator&&(simd_mask const& other) const {
    return simd_mask(vec_and(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask operator!() const {
    return simd_mask(vec_nand(m_value, simd_mask(true).get()));
  }
};

SIMD_ALWAYS_INLINE inline bool all_of(simd_mask<double, simd_abi::vsx> const& a) {
  auto const true_value = simd_mask<double, simd_abi::vsx>(true).get();
  return vec_all_eq(a.get(), true_value);
}

SIMD_ALWAYS_INLINE inline bool any_of(simd_mask<double, simd_abi::vsx> const& a) {
  auto const true_value = simd_mask<double, simd_abi::vsx>(true).get();
  return vec_any_eq(a.get(), true_value);
}

template <>
class simd<double, simd_abi::vsx> {
  __vector double m_value;
 public:
  using value_type = double;
  using abi_type = simd_abi::vsx;
  using mask_type = simd_mask<double, abi_type>;
  using storage_type = simd_storage<double, abi_type>;
  SIMD_ALWAYS_INLINE inline simd() = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 2; }
  SIMD_ALWAYS_INLINE inline simd(double value)
    :m_value(vec_splats(value))
  {}
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
  SIMD_ALWAYS_INLINE inline simd(double const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_ALWAYS_INLINE inline constexpr simd(__vector double const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline simd operator*(simd const& other) const {
    return simd(vec_mul(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator/(simd const& other) const {
    return simd(vec_div(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator+(simd const& other) const {
    return simd(vec_add(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator-(simd const& other) const {
    return simd(vec_sub(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd operator-() const {
    // return simd(vec_neg(m_value)); some GCC versions dont have this
    return simd(0.0) - (*this);
  }
  SIMD_ALWAYS_INLINE inline void copy_from(double const* ptr, element_aligned_tag) {
    m_value = vec_vsx_ld(0, ptr);
  }
  SIMD_ALWAYS_INLINE inline void copy_to(double* ptr, element_aligned_tag) const {
    vec_vsx_st(m_value, 0, ptr);
  }
  SIMD_ALWAYS_INLINE inline constexpr __vector double get() const { return m_value; }
  SIMD_ALWAYS_INLINE inline simd_mask<double, simd_abi::vsx> operator<(simd const& other) const {
    return simd_mask<double, simd_abi::vsx>(vec_cmplt(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask<double, simd_abi::vsx> operator==(simd const& other) const {
    return simd_mask<double, simd_abi::vsx>(vec_cmpeq(m_value, other.m_value));
  }
};

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::vsx> abs(simd<double, simd_abi::vsx> const& a) {
  return simd<double, simd_abi::vsx>(vec_abs(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::vsx> sqrt(simd<double, simd_abi::vsx> const& a) {
  return simd<double, simd_abi::vsx>(vec_sqrt(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::vsx> fma(
    simd<double, simd_abi::vsx> const& a,
    simd<double, simd_abi::vsx> const& b,
    simd<double, simd_abi::vsx> const& c) {
  return simd<double, simd_abi::vsx>(vec_madd(a.get(), b.get(), c.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::vsx> max(
    simd<double, simd_abi::vsx> const& a, simd<double, simd_abi::vsx> const& b) {
  return simd<double, simd_abi::vsx>(vec_max(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::vsx> min(
    simd<double, simd_abi::vsx> const& a, simd<double, simd_abi::vsx> const& b) {
  return simd<double, simd_abi::vsx>(vec_min(a.get(), b.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::vsx> choose(
    simd_mask<double, simd_abi::vsx> const& a, simd<double, simd_abi::vsx> const& b, simd<double, simd_abi::vsx> const& c) {
  return simd<double, simd_abi::vsx>(vec_sel(c.get(), b.get(), a.get()));
}

#endif

#ifdef __CUDACC__

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
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE
  simd_mask(bool value)
    :m_value(value)
  {}
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE constexpr
  bool get() const {
    return m_value;
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE 
  simd_mask operator||(simd_mask const& other) const {
    return m_value || other.m_value;
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE
  simd_mask operator&&(simd_mask const& other) const {
    return m_value && other.m_value;
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE
  simd_mask operator!() const {
    return !m_value;
  }
};

template <class T, int N>
SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE
bool all_of(simd_mask<T, simd_abi::cuda_warp<N>> const& a) {
  return bool(__all_sync(simd_abi::cuda_warp<N>::mask(), int(a.get())));
}

template <class T, int N>
SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE
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
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE simd(T value)
    :m_value(value)
  {}
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE 
  simd(storage_type const& value) {
    copy_from(value.data(), element_aligned_tag());
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE
  simd& operator=(storage_type const& value) {
    copy_from(value.data(), element_aligned_tag());
    return *this;
  }
  template <class Flags>
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE simd(T const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE simd operator*(simd const& other) const {
    return simd(m_value * other.m_value);
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE simd operator/(simd const& other) const {
    return simd(m_value / other.m_value);
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE simd operator+(simd const& other) const {
    return simd(m_value + other.m_value);
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE simd operator-(simd const& other) const {
    return simd(m_value - other.m_value);
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE simd operator-() const {
    return simd(-m_value);
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE void copy_from(T const* ptr, element_aligned_tag) {
    m_value = ptr[threadIdx.x];
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE void copy_to(T* ptr, element_aligned_tag) const {
    ptr[threadIdx.x] = m_value;
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE T get() const {
    return m_value;
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE
  mask_type operator<(simd const& other) const {
    return mask_type(m_value < other.m_value);
  }
  SIMD_CUDA_ALWAYS_INLINE SIMD_DEVICE
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

#endif

template <class T>
class simd_size {
  public:
  static constexpr int value = 1;
};

template <class T, class Abi>
class simd_size<simd<T, Abi>> {
  public:
  static constexpr int value = simd<T, Abi>::size();
};

namespace simd_abi {

#if defined(__CUDACC__)
using native = scalar;
#elif defined(__AVX512F__)
using native = avx512;
#elif defined(__AVX__)
using native = avx;
#elif defined(__SSE2__)
using native = sse;
#elif defined(__ARM_NEON)
using native = neon;
#elif defined(__VSX__)
using native = vsx;
#else
using native = pack<8>;
#endif

}

template <class T>
using native_simd = simd<T, simd_abi::native>;

}

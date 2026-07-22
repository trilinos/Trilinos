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

#include <cmath>
#include <cstdint>

#ifndef SIMD_ALWAYS_INLINE
#define SIMD_ALWAYS_INLINE [[gnu::always_inline]]
#endif

#if defined( __CUDACC__ )
#define SIMD_CUDA_ALWAYS_INLINE __forceinline__
#endif

#if defined( __HIPCC__ )
#define SIMD_HIP_ALWAYS_INLINE __forceinline__
#endif


#if defined( __CUDACC__) || defined( __HIPCC__ )
#define SIMD_HOST_DEVICE __host__ __device__
#else
#define SIMD_HOST_DEVICE
#endif

#if defined (__CUDACC__) || defined( __HIPCC__ )
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

SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline double copysign(double a, double b) {
  return std::copysign(a, b);
}

SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline float copysign(float a, float b) {
  return std::copysignf(a, b);
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi> multiplysign(simd<T, Abi> a, simd<T, Abi> b) {
  T tmp_a[simd<T, Abi>::size()];
  T tmp_b[simd<T, Abi>::size()];
  a.copy_to(tmp_a, element_aligned_tag());
  b.copy_to(tmp_b, element_aligned_tag());
  for (int i = 0; i < simd<T, Abi>::size(); ++i) tmp_a[i] = tmp_a[i]*::SIMD_NAMESPACE::copysign(static_cast<T>(1.0), tmp_b[i]);
  a.copy_from(tmp_a, element_aligned_tag());
  return a;
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi> copysign(simd<T, Abi> a, simd<T, Abi> b) {
  T tmp_a[simd<T, Abi>::size()];
  T tmp_b[simd<T, Abi>::size()];
  a.copy_to(tmp_a, element_aligned_tag());
  b.copy_to(tmp_b, element_aligned_tag());
  for (int i = 0; i < simd<T, Abi>::size(); ++i) tmp_a[i] = ::SIMD_NAMESPACE::copysign(tmp_a[i], tmp_b[i]);
  a.copy_from(tmp_a, element_aligned_tag());
  return a;
}

template <class T, class Abi>
SIMD_ALWAYS_INLINE SIMD_HOST_DEVICE inline simd<T, Abi> abs(simd<T, Abi> a) {
  T tmp[simd<T, Abi>::size()];
  a.copy_to(tmp, element_aligned_tag());
  for (int i = 0; i < simd<T, Abi>::size(); ++i) tmp[i] = std::abs(tmp[i]);
  a.copy_from(tmp, element_aligned_tag());
  return a;
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
  SIMD_ALWAYS_INLINE explicit inline
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

}

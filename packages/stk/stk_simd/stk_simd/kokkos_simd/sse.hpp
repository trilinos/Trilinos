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

#ifdef __SSE__
#include <xmmintrin.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#if defined(__FMA__) || defined(__AVX2__)
#include <immintrin.h>
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

namespace SIMD_NAMESPACE {

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
  SIMD_ALWAYS_INLINE inline simd(
      float a, float b, float c, float d)
    :m_value(_mm_setr_ps(a, b, c, d))
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
  SIMD_ALWAYS_INLINE inline simd(float const* ptr, Flags /*flags*/)
    :m_value(_mm_loadu_ps(ptr))
  {}
  SIMD_ALWAYS_INLINE inline simd(float const* ptr, int stride)
    :simd(ptr[0], ptr[stride], ptr[2*stride], ptr[3*stride])
  {}
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

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::sse> multiplysign(simd<float, simd_abi::sse> const& a, simd<float, simd_abi::sse> const& b) {
  __m128 const sign_mask = _mm_set1_ps(-0.);
  return simd<float, simd_abi::sse>(_mm_xor_ps(a.get(), _mm_and_ps(sign_mask, b.get())));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::sse> copysign(simd<float, simd_abi::sse> const& a, simd<float, simd_abi::sse> const& b) {
  __m128 const sign_mask = _mm_set1_ps(-0.);
  return simd<float, simd_abi::sse>(_mm_xor_ps(_mm_andnot_ps(sign_mask, a.get()), _mm_and_ps(sign_mask, b.get())));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::sse> abs(simd<float, simd_abi::sse> const& a) {
  __m128 const sign_mask = _mm_set1_ps(-0.f);  // -0.f = 1 << 31
  return simd<float, simd_abi::sse>(_mm_andnot_ps(sign_mask, a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::sse> sqrt(simd<float, simd_abi::sse> const& a) {
  return simd<float, simd_abi::sse>(_mm_sqrt_ps(a.get()));
}

#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
SIMD_ALWAYS_INLINE inline simd<float, simd_abi::sse> cbrt(simd<float, simd_abi::sse> const& a) {
  return simd<float, simd_abi::sse>(_mm_cbrt_ps(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::sse> exp(simd<float, simd_abi::sse> const& a) {
  return simd<float, simd_abi::sse>(_mm_exp_ps(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::sse> log(simd<float, simd_abi::sse> const& a) {
  return simd<float, simd_abi::sse>(_mm_log_ps(a.get()));
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
  SIMD_ALWAYS_INLINE inline simd(simd const&) = default;
  SIMD_ALWAYS_INLINE inline simd(simd&&) = default;
  SIMD_ALWAYS_INLINE inline simd& operator=(simd const&) = default;
  SIMD_ALWAYS_INLINE inline simd& operator=(simd&&) = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 2; }
  SIMD_ALWAYS_INLINE inline simd(double value)
    :m_value(_mm_set1_pd(value))
  {}
  SIMD_ALWAYS_INLINE inline simd(double a, double b)
    :m_value(_mm_setr_pd(a, b))
  {}
  SIMD_ALWAYS_INLINE inline
  simd(storage_type const& value) {
    copy_from(value.data(), element_aligned_tag());
  }
#ifdef STK_VOLATILE_SIMD
  SIMD_ALWAYS_INLINE inline
  simd(simd const volatile& value)
    :m_value(value.m_value)
  {}
#endif
  SIMD_ALWAYS_INLINE inline
  simd& operator=(storage_type const& value) {
    copy_from(value.data(), element_aligned_tag());
    return *this;
  }
  template <class Flags>
  SIMD_ALWAYS_INLINE inline simd(double const* ptr, Flags /*flags*/)
    :m_value(_mm_loadu_pd(ptr))
  {}
  SIMD_ALWAYS_INLINE inline simd(double const* ptr, int stride)
    :simd(ptr[0], ptr[stride])
  {}
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
#ifdef STK_VOLATILE_SIMD
  SIMD_ALWAYS_INLINE inline void plus_equals(simd const volatile& other) volatile {
    m_value = _mm_add_pd(m_value, other.m_value);
  }
#endif
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

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::sse> multiplysign(simd<double, simd_abi::sse> const& a, simd<double, simd_abi::sse> const& b) {
  __m128d const sign_mask = _mm_set1_pd(-0.);
  return simd<double, simd_abi::sse>(_mm_xor_pd(a.get(), _mm_and_pd(sign_mask, b.get())));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::sse> copysign(simd<double, simd_abi::sse> const& a, simd<double, simd_abi::sse> const& b) {
  __m128d const sign_mask = _mm_set1_pd(-0.);
  return simd<double, simd_abi::sse>(_mm_xor_pd(_mm_andnot_pd(sign_mask, a.get()), _mm_and_pd(sign_mask, b.get())));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::sse> abs(simd<double, simd_abi::sse> const& a) {
  __m128d const sign_mask = _mm_set1_pd(-0.);  // -0. = 1 << 63
  return simd<double, simd_abi::sse>(_mm_andnot_pd(sign_mask, a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::sse> sqrt(simd<double, simd_abi::sse> const& a) {
  return simd<double, simd_abi::sse>(_mm_sqrt_pd(a.get()));
}

#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
SIMD_ALWAYS_INLINE inline simd<double, simd_abi::sse> cbrt(simd<double, simd_abi::sse> const& a) {
  return simd<double, simd_abi::sse>(_mm_cbrt_pd(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::sse> exp(simd<double, simd_abi::sse> const& a) {
  return simd<double, simd_abi::sse>(_mm_exp_pd(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::sse> log(simd<double, simd_abi::sse> const& a) {
  return simd<double, simd_abi::sse>(_mm_log_pd(a.get()));
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

}

#endif

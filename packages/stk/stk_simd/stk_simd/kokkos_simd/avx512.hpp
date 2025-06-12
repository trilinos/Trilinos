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

#ifdef __AVX512F__

#include <immintrin.h>

namespace SIMD_NAMESPACE {

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
  static const __mmask16 false_value(-std::int16_t(false));
  return _kortestc_mask16_u8(a.get(), false_value);
}

SIMD_ALWAYS_INLINE inline bool any_of(simd_mask<float, simd_abi::avx512> const& a) {
  static const __mmask16 false_value(-std::int16_t(false));
  return !_kortestc_mask16_u8(~a.get(), false_value);
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
  SIMD_ALWAYS_INLINE inline simd(
      float a, float b, float c, float d,
      float e, float f, float g, float h,
      float i, float j, float k, float l,
      float m, float n, float o, float p)
    :m_value(_mm512_setr_ps(
          a, b, c, d, e, f, g, h,
          i, j, k, l, m, n, o, p))
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
    :m_value(_mm512_loadu_ps(ptr))
  {}
  SIMD_ALWAYS_INLINE inline simd(float const* ptr, int stride)
    :simd(ptr[0],        ptr[stride],   ptr[2*stride], ptr[3*stride],
          ptr[4*stride], ptr[5*stride], ptr[6*stride], ptr[7*stride],
          ptr[8*stride], ptr[9*stride], ptr[10*stride], ptr[11*stride],
          ptr[12*stride], ptr[13*stride], ptr[14*stride], ptr[15*stride])
  {}
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
    return simd_mask<float, simd_abi::avx512>(_mm512_cmplt_ps_mask(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask<float, simd_abi::avx512> operator==(simd const& other) const {
    return simd_mask<float, simd_abi::avx512>(_mm512_cmpeq_ps_mask(m_value, other.m_value));
  }
};

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx512> multiplysign(simd<float, simd_abi::avx512> const& a, simd<float, simd_abi::avx512> const& b) {
  static const __m512i sign_mask = reinterpret_cast<__m512i>(simd<float, simd_abi::avx512>(-0.0).get());
  return simd<float, simd_abi::avx512>(
      reinterpret_cast<__m512>(_mm512_xor_epi32(
          reinterpret_cast<__m512i>(a.get()), 
          _mm512_and_epi32(sign_mask, reinterpret_cast<__m512i>(b.get()))
          ))
      );
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx512> copysign(simd<float, simd_abi::avx512> const& a, simd<float, simd_abi::avx512> const& b) {
  static const __m512i sign_mask = reinterpret_cast<__m512i>(simd<float, simd_abi::avx512>(-0.0).get());
  return simd<float, simd_abi::avx512>(
      reinterpret_cast<__m512>(_mm512_xor_epi32(
          _mm512_andnot_epi32(sign_mask, reinterpret_cast<__m512i>(a.get())),
          _mm512_and_epi32(sign_mask, reinterpret_cast<__m512i>(b.get()))
          ))
      );
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx512> abs(simd<float, simd_abi::avx512> const& a) {
  __m512 const rhs = a.get();
  return reinterpret_cast<__m512>(_mm512_and_epi32(reinterpret_cast<__m512i>(rhs), _mm512_set1_epi32(0x7fffffff)));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx512> sqrt(simd<float, simd_abi::avx512> const& a) {
  return simd<float, simd_abi::avx512>(_mm512_sqrt_ps(a.get()));
}

#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx512> cbrt(simd<float, simd_abi::avx512> const& a) {
  return simd<float, simd_abi::avx512>(_mm512_cbrt_ps(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx512> exp(simd<float, simd_abi::avx512> const& a) {
  return simd<float, simd_abi::avx512>(_mm512_exp_ps(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx512> log(simd<float, simd_abi::avx512> const& a) {
  return simd<float, simd_abi::avx512>(_mm512_log_ps(a.get()));
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
    return simd_mask(static_cast<__mmask8>(_mm512_kor(m_value, other.m_value)));
  }
  SIMD_ALWAYS_INLINE simd_mask operator&&(simd_mask const& other) const {
    return simd_mask(static_cast<__mmask8>(_mm512_kand(m_value, other.m_value)));
  }
  SIMD_ALWAYS_INLINE simd_mask operator!() const {
    static const __mmask8 true_value(simd_mask<double, simd_abi::avx512>(true).get());
    return simd_mask(static_cast<__mmask8>(_mm512_kxor(true_value, m_value)));
  }
};

SIMD_ALWAYS_INLINE inline bool all_of(simd_mask<double, simd_abi::avx512> const& a) {
  static const __mmask16 false_value(-std::int16_t(false));
  const __mmask16 a_value(0xFF00 | a.get());
  return _kortestc_mask16_u8(a_value, false_value);
}

SIMD_ALWAYS_INLINE inline bool any_of(simd_mask<double, simd_abi::avx512> const& a) {
  static const __mmask16 false_value(-std::int16_t(false));
  const __mmask16 a_value(0x0000 | a.get());
  return !_kortestc_mask16_u8(~a_value, false_value);
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
  SIMD_ALWAYS_INLINE inline simd(simd const&) = default;
  SIMD_ALWAYS_INLINE inline simd(simd&&) = default;
  SIMD_ALWAYS_INLINE inline simd& operator=(simd const&) = default;
  SIMD_ALWAYS_INLINE inline simd& operator=(simd&&) = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 8; }
  SIMD_ALWAYS_INLINE inline simd(double value)
    :m_value(_mm512_set1_pd(value))
  {}
  SIMD_ALWAYS_INLINE inline simd(
      double a, double b, double c, double d,
      double e, double f, double g, double h)
    :m_value(_mm512_setr_pd(a, b, c, d, e, f, g, h))
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
    :m_value(_mm512_loadu_pd(ptr))
  {}
  SIMD_ALWAYS_INLINE inline simd(double const* ptr, int stride)
    :simd(ptr[0],        ptr[stride],   ptr[2*stride], ptr[3*stride],
          ptr[4*stride], ptr[5*stride], ptr[6*stride], ptr[7*stride])
  {}
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
#ifdef STK_VOLATILE_SIMD
  SIMD_ALWAYS_INLINE inline void plus_equals(simd const volatile& other) volatile {
    m_value = _mm512_add_pd(m_value, other.m_value);
  }
#endif
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
    return simd_mask<double, simd_abi::avx512>(_mm512_cmplt_pd_mask(m_value, other.m_value));
  }
  SIMD_ALWAYS_INLINE inline simd_mask<double, simd_abi::avx512> operator==(simd const& other) const {
    return simd_mask<double, simd_abi::avx512>(_mm512_cmpeq_pd_mask(m_value, other.m_value));
  }
};

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx512> multiplysign(simd<double, simd_abi::avx512> const& a, simd<double, simd_abi::avx512> const& b) {
  static const __m512i sign_mask = reinterpret_cast<__m512i>(simd<double, simd_abi::avx512>(-0.0).get());
  return simd<double, simd_abi::avx512>(
      reinterpret_cast<__m512d>(_mm512_xor_epi64(
          reinterpret_cast<__m512i>(a.get()), 
          _mm512_and_epi64(sign_mask, reinterpret_cast<__m512i>(b.get()))
          ))
      );
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx512> copysign(simd<double, simd_abi::avx512> const& a, simd<double, simd_abi::avx512> const& b) {
  static const __m512i sign_mask = reinterpret_cast<__m512i>(simd<double, simd_abi::avx512>(-0.0).get());
  return simd<double, simd_abi::avx512>(
      reinterpret_cast<__m512d>(_mm512_xor_epi64(
          _mm512_andnot_epi64(sign_mask, reinterpret_cast<__m512i>(a.get())),
          _mm512_and_epi64(sign_mask, reinterpret_cast<__m512i>(b.get()))
          ))
      );
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx512> abs(simd<double, simd_abi::avx512> const& a) {
  __m512d const rhs = a.get();
  return reinterpret_cast<__m512d>(_mm512_and_epi64(_mm512_set1_epi64(0x7FFFFFFFFFFFFFFF),
        reinterpret_cast<__m512i>(rhs)));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx512> sqrt(simd<double, simd_abi::avx512> const& a) {
  return simd<double, simd_abi::avx512>(_mm512_sqrt_pd(a.get()));
}

#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx512> cbrt(simd<double, simd_abi::avx512> const& a) {
  return simd<double, simd_abi::avx512>(_mm512_cbrt_pd(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx512> exp(simd<double, simd_abi::avx512> const& a) {
  return simd<double, simd_abi::avx512>(_mm512_exp_pd(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx512> log(simd<double, simd_abi::avx512> const& a) {
  return simd<double, simd_abi::avx512>(_mm512_log_pd(a.get()));
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

}

#endif




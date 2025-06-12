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

#ifdef __AVX__

#include <immintrin.h>

namespace SIMD_NAMESPACE {

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
  SIMD_ALWAYS_INLINE inline simd(
      float a, float b, float c, float d,
      float e, float f, float g, float h)
    :m_value(_mm256_setr_ps(a, b, c, d, e, f, g, h))
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
    :m_value(_mm256_loadu_ps(ptr))
  {}
  SIMD_ALWAYS_INLINE inline simd(float const* ptr, int stride)
    :simd(ptr[0],        ptr[stride],   ptr[2*stride], ptr[3*stride],
          ptr[4*stride], ptr[5*stride], ptr[6*stride], ptr[7*stride])
  {}
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

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx> multiplysign(simd<float, simd_abi::avx> const& a, simd<float, simd_abi::avx> const& b) {
  __m256 const sign_mask = _mm256_set1_ps(-0.f);
  return simd<float, simd_abi::avx>(_mm256_xor_ps(a.get(), _mm256_and_ps(sign_mask, b.get())));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx> copysign(simd<float, simd_abi::avx> const& a, simd<float, simd_abi::avx> const& b) {
  __m256 const sign_mask = _mm256_set1_ps(-0.);
  return simd<float, simd_abi::avx>(_mm256_xor_ps(_mm256_andnot_ps(sign_mask, a.get()), _mm256_and_ps(sign_mask, b.get())));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx> abs(simd<float, simd_abi::avx> const& a) {
  __m256 sign_mask = _mm256_set1_ps(-0.f);  // -0.f = 1 << 31
  return simd<float, simd_abi::avx>(_mm256_andnot_ps(sign_mask, a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx> sqrt(simd<float, simd_abi::avx> const& a) {
  return simd<float, simd_abi::avx>(_mm256_sqrt_ps(a.get()));
}

#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx> cbrt(simd<float, simd_abi::avx> const& a) {
  return simd<float, simd_abi::avx>(_mm256_cbrt_ps(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx> exp(simd<float, simd_abi::avx> const& a) {
  return simd<float, simd_abi::avx>(_mm256_exp_ps(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<float, simd_abi::avx> log(simd<float, simd_abi::avx> const& a) {
  return simd<float, simd_abi::avx>(_mm256_log_ps(a.get()));
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
  SIMD_ALWAYS_INLINE inline simd(simd const&) = default;
  SIMD_ALWAYS_INLINE inline simd(simd&&) = default;
  SIMD_ALWAYS_INLINE inline simd& operator=(simd const&) = default;
  SIMD_ALWAYS_INLINE inline simd& operator=(simd&&) = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 4; }
  SIMD_ALWAYS_INLINE inline simd(double value)
    :m_value(_mm256_set1_pd(value))
  {}
  SIMD_ALWAYS_INLINE inline simd(
      double a, double b, double c, double d)
    :m_value(_mm256_setr_pd(a, b, c, d))
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
  SIMD_ALWAYS_INLINE inline simd(double const* ptr, Flags flags)
    :m_value(_mm256_loadu_pd(ptr))
  {}
  SIMD_ALWAYS_INLINE inline simd(double const* ptr, int stride)
    :simd(ptr[0], ptr[stride], ptr[2*stride], ptr[3*stride])
  {}
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
#ifdef STK_VOLATILE_SIMD
  SIMD_ALWAYS_INLINE inline void plus_equals(simd const volatile& other) volatile {
    m_value = _mm256_add_pd(m_value, other.m_value);
  }
#endif
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

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx> multiplysign(simd<double, simd_abi::avx> const& a, simd<double, simd_abi::avx> const& b) {
  __m256d const sign_mask = _mm256_set1_pd(-0.f);
  return simd<double, simd_abi::avx>(_mm256_xor_pd(a.get(), _mm256_and_pd(sign_mask, b.get())));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx> copysign(simd<double, simd_abi::avx> const& a, simd<double, simd_abi::avx> const& b) {
  __m256d const sign_mask = _mm256_set1_pd(-0.f);
  return simd<double, simd_abi::avx>(_mm256_xor_pd(_mm256_andnot_pd(sign_mask, a.get()), _mm256_and_pd(sign_mask, b.get())));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx> abs(simd<double, simd_abi::avx> const& a) {
  __m256d const sign_mask = _mm256_set1_pd(-0.f);  // -0.f = 1 << 31
  return simd<double, simd_abi::avx>(_mm256_andnot_pd(sign_mask, a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx> sqrt(simd<double, simd_abi::avx> const& a) {
  return simd<double, simd_abi::avx>(_mm256_sqrt_pd(a.get()));
}

#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx> cbrt(simd<double, simd_abi::avx> const& a) {
  return simd<double, simd_abi::avx>(_mm256_cbrt_pd(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx> exp(simd<double, simd_abi::avx> const& a) {
  return simd<double, simd_abi::avx>(_mm256_exp_pd(a.get()));
}

SIMD_ALWAYS_INLINE inline simd<double, simd_abi::avx> log(simd<double, simd_abi::avx> const& a) {
  return simd<double, simd_abi::avx>(_mm256_log_pd(a.get()));
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

}

#endif

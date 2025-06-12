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

#ifdef __ARM_NEON

#include <arm_neon.h>

namespace SIMD_NAMESPACE {

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
  SIMD_ALWAYS_INLINE inline simd(float a, float b, float c, float d)
    :m_value((float32x4_t){a, b, c, d})
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
  SIMD_ALWAYS_INLINE inline simd(float const* ptr, int stride)
    :simd(ptr[0], ptr[stride], ptr[2*stride], ptr[3*stride])
  {}
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
  SIMD_ALWAYS_INLINE inline simd(simd const&) = default;
  SIMD_ALWAYS_INLINE inline simd(simd&&) = default;
  SIMD_ALWAYS_INLINE inline simd& operator=(simd const&) = default;
  SIMD_ALWAYS_INLINE inline simd& operator=(simd&&) = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 2; }
  SIMD_ALWAYS_INLINE inline simd(double value)
    :m_value(vdupq_n_f64(value))
  {}
  SIMD_ALWAYS_INLINE inline simd(double a, double b)
    :m_value((float64x2_t){a, b})
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
  SIMD_ALWAYS_INLINE inline simd(double const* ptr, Flags flags) {
    copy_from(ptr, flags);
  }
  SIMD_ALWAYS_INLINE inline simd(double const* ptr, int stride)
    :simd(ptr[0], ptr[stride])
  {}
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
#ifdef STK_VOLATILE_SIMD
  SIMD_ALWAYS_INLINE inline void plus_equals(simd const volatile& other) volatile {
    m_value = vaddq_f64(m_value, other.m_value);
  }
#endif
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

}

#endif

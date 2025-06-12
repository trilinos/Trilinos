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

#ifdef __VSX__
#include <altivec.h>
// undefine the really dangerous macros from this file
#undef vector
#undef pixel
#undef bool
#endif


#if defined(__VSX__) && (!defined(__CUDACC__))

namespace SIMD_NAMESPACE {

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
  SIMD_ALWAYS_INLINE inline simd(float a, float b, float c, float d)
    :m_value((__vector float){a, b, c, d})
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
  __vector __bool long long m_value;
  using ll_t = long long;
  using ull_t = unsigned long long;
 public:
  using value_type = bool;
  using simd_type = simd_mask<double, simd_abi::vsx>;
  using abi_type = simd_abi::vsx;
  SIMD_ALWAYS_INLINE inline simd_mask() = default;
  SIMD_ALWAYS_INLINE inline simd_mask(bool value)
    :m_value{ull_t(-ll_t(value)), ull_t(-ll_t(value))}
  {}
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 2; }
  SIMD_ALWAYS_INLINE inline constexpr simd_mask(__vector __bool long long const& value_in)
    :m_value(value_in)
  {}
  SIMD_ALWAYS_INLINE inline constexpr __vector __bool long long get() const { return m_value; }
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
  SIMD_ALWAYS_INLINE inline simd(simd const&) = default;
  SIMD_ALWAYS_INLINE inline simd(simd&&) = default;
  SIMD_ALWAYS_INLINE inline simd& operator=(simd const&) = default;
  SIMD_ALWAYS_INLINE inline simd& operator=(simd&&) = default;
  SIMD_ALWAYS_INLINE inline static constexpr int size() { return 2; }
  SIMD_ALWAYS_INLINE inline simd(double value)
    :m_value(vec_splats(value))
  {}
  SIMD_ALWAYS_INLINE inline simd(double a, double b)
    :m_value((__vector double){a, b})
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
#ifdef STK_VOLATILE_SIMD
  SIMD_ALWAYS_INLINE inline void plus_equals(simd const volatile& other) volatile {
    m_value = vec_add(m_value, other.m_value);
  }
#endif
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

}

#endif

// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// IWYU pragma: private; include "Simd.hpp"

#ifndef STK_INCLUDE_ONLY_STK_SIMD_HEADER
static_assert(false, "Do not include simd impl files directly. Only include stk_simd/Simd.hpp");
#endif

#include <type_traits>

#ifndef STK_SIMD_AVX512SPECIALIZATIONS_HPP
#define STK_SIMD_AVX512SPECIALIZATIONS_HPP

#if defined(__AVX512F__) && !defined(USE_STK_SIMD_NONE)
#define STK_USE_AVX512_SPECIALIZATION
#endif

namespace stk::simd::impl::specialized
{

#ifdef STK_USE_AVX512_SPECIALIZATION

template <class T, class simd_type = typename SimdT<T>::type, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE constexpr auto avx512_make_mask(const int numValid)
{
  using mask_type = std::conditional_t<std::is_same_v<Double, simd_type>, uint8_t, uint16_t>;
  assert(numValid >= 0 && numValid <= Traits<simd_type>::length);
  return static_cast<mask_type>((1 << numValid) - 1);
}

template <class T, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE auto avx512_make_index(const int offset)
{
  // offset is typically your "size" in your load_array() function.
  // The computation 7LL * offset etc. will be constant-folded when offset is compile-time constant.
  if constexpr (std::is_same_v<T, double>) {
    return _mm256_set_epi32(7 * offset, 6 * offset, 5 * offset, 4 * offset, 3 * offset, 2 * offset, offset, 0);
  } else {
    return _mm512_set_epi32(15 * offset, 14 * offset, 13 * offset, 12 * offset, 11 * offset, 10 * offset, 9 * offset,
        8 * offset, 7 * offset, 6 * offset, 5 * offset, 4 * offset, 3 * offset, 2 * offset, offset, 0);
  }
}

template <class T, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE auto avx512_packed_zero()
{
  if constexpr (std::is_same_v<T, double>) {
    return _mm512_setzero_pd();
  } else {
    return _mm512_setzero_ps();
  }
}

// FLOAT
STK_MATH_FORCE_INLINE Float avx512_masked_load(const __mmask16& mask, const float* const from)
{
  return Float(_mm512_maskz_loadu_ps(mask, from));
}

STK_MATH_FORCE_INLINE void avx512_masked_store(float* const dst, const __mmask16& mask, const Float& src)
{
  _mm512_mask_storeu_ps(dst, mask, src._data.get());
}

STK_MATH_FORCE_INLINE Float avx512_gather(const __m512i vindex, const float* const from)
{
  return Float(_mm512_i32gather_ps(vindex, from, sizeof(float)));
}

STK_MATH_FORCE_INLINE Float avx512_masked_gather(const __mmask16& mask, const __m512i& vindex, const float* const from)
{
  const auto zero = avx512_packed_zero<float>();
  return Float(_mm512_mask_i32gather_ps(zero, mask, vindex, from, sizeof(float)));
}

STK_MATH_FORCE_INLINE void avx512_scatter(float* const dst, const __m512i& vindex, const Float& src)
{
  _mm512_i32scatter_ps(dst, vindex, src._data.get(), sizeof(float));
}

STK_MATH_FORCE_INLINE void avx512_masked_scatter(
    float* const dst, const __mmask16& mask, const __m512i& vindex, const Float& src)
{
  _mm512_mask_i32scatter_ps(dst, mask, vindex, src._data.get(), sizeof(float));
}

// DOUBLE
STK_MATH_FORCE_INLINE Double avx512_masked_load(const __mmask8& mask, const double* const from)
{
  return Double(_mm512_maskz_loadu_pd(mask, from));
}

STK_MATH_FORCE_INLINE void avx512_masked_store(double* const dst, const __mmask8& mask, const Double& src)
{
  _mm512_mask_storeu_pd(dst, mask, src._data.get());
}

STK_MATH_FORCE_INLINE Double avx512_gather(const __m256i vindex, const double* const from)
{
  return Double(_mm512_i32gather_pd(vindex, from, sizeof(double)));
}

STK_MATH_FORCE_INLINE Double avx512_masked_gather(const __mmask8& mask, const __m256i& vindex, const double* const from)
{
  const auto zero = avx512_packed_zero<double>();
  return Double(_mm512_mask_i32gather_pd(zero, mask, vindex, from, sizeof(double)));
}

STK_MATH_FORCE_INLINE void avx512_scatter(double* const dst, const __m256i& vindex, const Double& src)
{
  _mm512_i32scatter_pd(dst, vindex, src._data.get(), sizeof(double));
}

STK_MATH_FORCE_INLINE void avx512_masked_scatter(
    double* const dst, const __mmask8& mask, const __m256i& vindex, const Double& src)
{
  _mm512_mask_i32scatter_pd(dst, mask, vindex, src._data.get(), sizeof(double));
}

#endif  // ifdef STK_USE_AVX512_SPECIALIZATION

}  // namespace stk::simd::impl::specialized

#endif

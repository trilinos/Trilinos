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

namespace stk::simd::impl
{

template <class T, class simd_type = typename SimdT<T>::type, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE simd_type load_part(const T* x, const int numValid)
{
  assert(numValid <= Traits<simd_type>::length);
  simd_type temp(0.0);
  for (int n = 0; n < numValid; ++n) {
    temp[n] = x[n];
  }
  return temp;
}

template <class T, class simd_type = typename SimdT<T>::type, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE simd_type load_part(const T* x, const int offset, const int numValid)
{
  assert(numValid <= Traits<simd_type>::length);
  simd_type temp(0.0);
  for (int n = 0; n < numValid; ++n) {
    temp[n] = x[n * offset];
  }
  return temp;
}

template <class T, class simd_type = typename SimdT<T>::type, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE void store_part(T* x, const simd_type& z, const int numValid)
{
  assert(numValid <= Traits<simd_type>::length);
  for (int n = 0; n < numValid; ++n) {
    x[n] = z[n];
  }
}

template <class T, class simd_type = typename SimdT<T>::type>
STK_MATH_FORCE_INLINE void store_part(T* x, const simd_type& z, const int offset, const int numValid)
{
  assert(numValid <= Traits<simd_type>::length);
  for (int n = 0; n < numValid; ++n) {
    x[n * offset] = z[n];
  }
}

template <int size, class T, class simd_type = typename SimdT<T>::type, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE void load_array(simd_type* const to, const T* const from)
{
  for (int i = 0; i < size; ++i) {
    to[i] = load(from + i, size);
  }
}

template <int size, class T, class simd_type = typename SimdT<T>::type, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE void load_array(simd_type* const to, const T* const from, const int numValid)
{
  for (int i = 0; i < size; ++i) {
    to[i] = impl::load_part(from + i, size, numValid);
  }
}

template <int size, class T, class simd_type = typename SimdT<T>::type, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE void store_array(T* const to, const simd_type* const from)
{
  for (int i = 0; i < size; ++i) {
    store(to + i, from[i], size);
  }
}

template <int size, class T, class simd_type = typename SimdT<T>::type, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE void store_array(T* const to, const simd_type* const from, const int numValid)
{
  for (int i = 0; i < size; ++i) {
    impl::store_part(to + i, from[i], size, numValid);
  }
}

namespace specialized
{

template <class T, class simd_type = typename SimdT<T>::type, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE void store_part(T* x, const simd_type& z, const int numValid)
{
#ifdef STK_USE_AVX512_SPECIALIZATION
  avx512_masked_store(x, avx512_make_mask<T>(numValid), z);
#else
  impl::store_part(x, z, numValid);
#endif
}

template <class T, class simd_type = typename SimdT<T>::type, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE void store_part(T* x, const simd_type& z, const int offset, const int numValid)
{
#ifdef STK_USE_AVX512_SPECIALIZATION
  if (numValid == Traits<simd_type>::length) {
    avx512_scatter(x, avx512_make_index<T>(offset), z);
  } else {
    avx512_masked_scatter(x, avx512_make_mask<T>(numValid), avx512_make_index<T>(offset), z);
  }
#else
  impl::store_part(x, z, offset, numValid);
#endif
}

template <class T, class simd_type = typename SimdT<T>::type, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE simd_type load_part(const T* x, const int numValid)
{
#ifdef STK_USE_AVX512_SPECIALIZATION
  return avx512_masked_load(avx512_make_mask<T>(numValid), x);
#else
  return impl::load_part(x, numValid);
#endif
}

template <class T, class simd_type = typename SimdT<T>::type, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE simd_type load_part(const T* x, const int offset, const int numValid)
{
#ifdef STK_USE_AVX512_SPECIALIZATION
  return avx512_masked_gather(avx512_make_mask<T>(numValid), avx512_make_index<T>(offset), x);
#else
  return impl::load_part(x, offset, numValid);
#endif
}

template <int size, class T, class simd_type = typename SimdT<T>::type, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE void load_array(simd_type* const to, const T* const from)
{
#ifdef STK_USE_AVX512_SPECIALIZATION
  auto vindex = avx512_make_index<T>(size);
  for (int i = 0; i < size; ++i) {
    to[i] = avx512_gather(vindex, from + i);
  }
#else
  impl::load_array<size>(to, from);
#endif
}

template <int size, class T, class simd_type = typename SimdT<T>::type, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE void store_array(T* const to, const simd_type* const from)
{
#ifdef STK_USE_AVX512_SPECIALIZATION
  auto vindex = avx512_make_index<T>(size);
  for (int i = 0; i < size; ++i) {
    avx512_scatter(to + i, vindex, from[i]);
  }
#else
  impl::store_array<size>(to, from);
#endif
}

template <int size, class T, class simd_type = typename SimdT<T>::type, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE void load_array(simd_type* const to, const T* const from, const int numValid)
{
#ifdef STK_USE_AVX512_SPECIALIZATION
  auto mask = avx512_make_mask<T>(numValid);
  auto vindex = avx512_make_index<T>(size);
  for (int i = 0; i < size; ++i) {
    to[i] = avx512_masked_gather(mask, vindex, from + i);
  }
#else
  impl::load_array<size>(to, from, numValid);
#endif
}

template <int size, class T, class simd_type = typename SimdT<T>::type, class = enable_if_supported_real_t<T>>
STK_MATH_FORCE_INLINE void store_array(T* const to, const simd_type* const from, const int numValid)
{
#ifdef STK_USE_AVX512_SPECIALIZATION
  auto mask = avx512_make_mask<T>(numValid);
  auto vindex = avx512_make_index<T>(size);
  for (int i = 0; i < size; ++i) {
    avx512_masked_scatter(to + i, mask, vindex, from[i]);
  }
#else
  impl::store_array<size>(to, from, numValid);
#endif
}

}  // namespace specialized

}  // namespace stk::simd::impl

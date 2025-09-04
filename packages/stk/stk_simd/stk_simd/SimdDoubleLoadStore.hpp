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

#ifndef STK_SIMD_DOUBLELOADSTORE_HPP
#define STK_SIMD_DOUBLELOADSTORE_HPP

namespace stk {
namespace simd
{
namespace impl {

template<int N>
STK_MATH_FORCE_INLINE
void store_strided(double* x, const simd::Double& z, const int offset)
{
  for(int i=0; i<N; ++i) x[i*offset] = z[i];
}

template<>
STK_MATH_FORCE_INLINE
void store_strided<2>(double* x, const simd::Double& z, const int offset)
{
  x[0] = z[0]; x[offset] = z[1];
}

template<>
STK_MATH_FORCE_INLINE
void store_strided<4>(double* x, const simd::Double& z, const int offset)
{
  x[0] = z[0]; x[offset] = z[1]; x[2*offset] = z[2]; x[3*offset] = z[3];
}

template<>
STK_MATH_FORCE_INLINE
void store_strided<8>(double* x, const simd::Double& z, const int offset)
{
#ifdef STK_USE_AVX512_SPECIALIZATION
  specialized::avx512_scatter(x, specialized::avx512_make_index<double>(offset), z);
#else
  x[0] = z[0]; x[offset] = z[1]; x[2*offset] = z[2]; x[3*offset] = z[3];
  x[4*offset] = z[4]; x[5*offset] = z[5]; x[6*offset] = z[6]; x[7*offset] = z[7];
#endif
}

} // namespace impl

STK_MATH_FORCE_INLINE simd::Double load_aligned(const double* x) {
  return simd::Double(SIMD_NAMESPACE::native_simd<double>(x, SIMD_NAMESPACE::element_aligned_tag()));
}

STK_MATH_FORCE_INLINE simd::Double load(const double* x) {
  //FIXME: supposed to be un-aligned load...
  return simd::Double(SIMD_NAMESPACE::native_simd<double>(x, SIMD_NAMESPACE::element_aligned_tag()));
}
    
STK_MATH_FORCE_INLINE simd::Double load(const double* x, const int offset) {
  return simd::Double(SIMD_NAMESPACE::native_simd<double>(x, offset));
}
  
STK_MATH_FORCE_INLINE void store_aligned(double* x, const simd::Double& z) {
  z._data.copy_to(x, SIMD_NAMESPACE::element_aligned_tag());
}

STK_MATH_FORCE_INLINE void store(double* x, const simd::Double& z) {
  //FIXME: supposed to be un-aligned store...
  z._data.copy_to(x, SIMD_NAMESPACE::element_aligned_tag());
}
  
STK_MATH_FORCE_INLINE void store(double* x, const simd::Double& z, const int offset) {
  impl::store_strided<ndoubles>(x, z, offset);
}

}  // namespace simd
} // namespace stk

#endif

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

#include "scalar.hpp"

#include "pack.hpp"

#if defined(__clang_major__) && (__clang_major__ >= 11)
#include "vector_size.hpp"
#endif

#ifndef SIMD_FORCE_SCALAR
#if defined( __CUDACC__ )
#include "cuda_warp.hpp"

#elif defined( __HIPCC__ )
#include "hip_wavefront.hpp"

#else

#ifdef __SSE__
#include "sse.hpp"
#endif

#ifdef __AVX__
#include "avx.hpp"
#endif

#ifdef __AVX512F__
#include "avx512.hpp"
#endif

#ifdef __ARM_NEON
#include "neon.hpp"
#endif

#ifdef __VSX__
#include "vsx.hpp"
#endif

#endif
#endif

namespace SIMD_NAMESPACE {

namespace simd_abi {

#if defined(SIMD_FORCE_SCALAR)
using native = scalar;
#elif defined(__CUDACC__)
using native = scalar;
#elif defined(__HIPCC__) 
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

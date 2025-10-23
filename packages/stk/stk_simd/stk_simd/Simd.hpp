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
// 

#ifndef STK_SIMD_SIMD_HPP
#define STK_SIMD_SIMD_HPP

#ifndef STK_INCLUDE_ONLY_STK_SIMD_HEADER
#define STK_INCLUDE_ONLY_STK_SIMD_HEADER
#endif

#include <stk_util/stk_config.h>

#include <cstdint>
#include <iostream>
#include <stk_math/StkMath.hpp>

#define STK_HAVE_SIMD
#define SIMD_NAMESPACE kokkos_simd

#ifndef SIMD_ALWAYS_INLINE
//currently necessary to avoid the 'always_inline' defined in simd.hpp
#define SIMD_ALWAYS_INLINE __attribute__((always_inline))
#endif

#ifdef USE_STK_SIMD_NONE
#define SIMD_FORCE_SCALAR
#endif

#include "kokkos_simd/simd.hpp"

namespace stk {
namespace simd {
inline constexpr int ndoubles = SIMD_NAMESPACE::native_simd<double>::size();
inline constexpr int nfloats = SIMD_NAMESPACE::native_simd<float>::size();
}
}

#include "SimdDouble.hpp"            // IWYU pragma: export
#include "SimdFloat.hpp"             // IWYU pragma: export
#include "SimdBool.hpp"              // IWYU pragma: export
#include "SimdBoolF.hpp"             // IWYU pragma: export
//
// has to be included after Double, Bool, Float, Boolf are defined
#include "Traits.hpp"  // IWYU pragma: export
//
#include "SimdAVX512Specializations.hpp"  // IWYU pragma: export
//
#include "SimdDoubleOperators.hpp"   // IWYU pragma: export
#include "SimdDoubleLoadStore.hpp"   // IWYU pragma: export
#include "SimdDoubleMath.hpp"        // IWYU pragma: export
//
#include "SimdFloatOperators.hpp"    // IWYU pragma: export
#include "SimdFloatLoadStore.hpp"    // IWYU pragma: export
#include "SimdFloatMath.hpp"         // IWYU pragma: export
//
#include "SimdLoadStoreImpl.hpp"  // IWYU pragma: export
//
#include <sys/resource.h>
#include <sys/time.h>

#include <Kokkos_Macros.hpp>

#include "stk_util/util/AlignedAllocator.hpp"

namespace stk {

inline double get_time_in_seconds() {
  timeval tp;
  struct timezone tz;
  gettimeofday(&tp, &tz);
  return (tp.tv_sec + tp.tv_usec/1000000.0);
}

namespace simd {

STK_MATH_FORCE_INLINE double reduce_sum(const Double& x) {
#ifdef STK_USE_AVX512_SPECIALIZATION
  return _mm512_reduce_add_pd(x._data.get());
#else
  double sum = x[0];
  for (int i=1; i<ndoubles; ++i) {
    sum += x[i];
  }
  return sum;
#endif
}

STK_MATH_FORCE_INLINE float reduce_sum(const Float& x) {
  double sum = x[0];
  for (int i=1; i<nfloats; ++i) {
    sum += x[i];
  }
  return sum;
}

STK_MATH_FORCE_INLINE double reduce_max(const Double& x) {
#ifdef STK_USE_AVX512_SPECIALIZATION
  return _mm512_reduce_max_pd(x._data.get());
#else
  double max_val = x[0];
  for (int i = 1; i < ndoubles; ++i) {
    if (x[i] > max_val) {
      max_val = x[i];
    }
  }
  return max_val;
#endif
}

STK_MATH_FORCE_INLINE float reduce_max(const Float& x) {
  float max = x[0];
  for (int i=1; i<nfloats; ++i) {
    max = max > x[i] ? max : x[i];
  }
  return max;
}

STK_MATH_FORCE_INLINE double reduce_min(const Double& x) {
#ifdef STK_USE_AVX512_SPECIALIZATION
  return _mm512_reduce_min_pd(x._data.get());
#else
  double min_val = x[0];
  for (int i = 1; i < ndoubles; ++i) {
    if (x[i] < min_val) {
      min_val = x[i];
    }
  }
  return min_val;
#endif
}

STK_MATH_FORCE_INLINE float reduce_min(const Float& x) {
  float max = x[0];
  for (int i=1; i<nfloats; ++i) {
    max = max < x[i] ? max : x[i];
  }
  return max;
}

//
//  These are helper functions for processing remainder loops.  The end pieces of SIMD loops
//  that do not have a full 'stk::ndoubles' entries.  'numValid' are the number of actual
//  entries to be processed.  Assumed that numValid is at least one and is less than
//  stk::ndoubles
//

// double versions

STK_MATH_FORCE_INLINE void store_part(double* x, const Double& z, const int numValid) {
  impl::specialized::store_part(x, z, numValid);
}

STK_MATH_FORCE_INLINE void store_part(double* x, const Double& z, const int offset, const int numValid) {
  impl::specialized::store_part(x, z, offset, numValid);
}

STK_MATH_FORCE_INLINE Double load_part(const double* x, const int numValid) {
  return impl::specialized::load_part(x, numValid);
}

STK_MATH_FORCE_INLINE Double load_part(const double* x, const int offset, const int numValid) {
  return impl::specialized::load_part(x, offset, numValid);
}

STK_MATH_FORCE_INLINE const double& get_data(const Double& z, int index) {
  assert(index < ndoubles);
  return z[index];
}

STK_MATH_FORCE_INLINE double& get_data(Double& z, int index) {
  assert(index < ndoubles);
  return z[index];
}

STK_MATH_FORCE_INLINE void set_data(Double& z, int index, const double val) {
  assert(index < ndoubles);
  z[index] = val;
}

inline std::ostream& operator << (std::ostream& output, const Double& a) {
  output << "[ ";
  for ( int i=0; i<ndoubles; ++i ) {
    output << a[i];
    if ( i < ndoubles - 1) {
      output << ", ";
    }
  }
  output << " ]";
  return output;
}

STK_MATH_FORCE_INLINE double reduce_sum(const Double& x, const int sumNum) {
  double ret = x[0];
  for (int i=1; i < sumNum; ++i) {
    ret += x[i];
  }
  return ret;
}

STK_MATH_FORCE_INLINE double reduce_max(const Double& x, const int sumNum) {
  double max = x[0];
  for (int i=1; i<sumNum; ++i) {
    max = max > x[i] ? max : x[i];
  }
  return max;
}

STK_MATH_FORCE_INLINE double reduce_min(const Double& x, const int sumNum) {
  double min = x[0];
  for (int i=1; i<sumNum; ++i) {
    min = min < x[i] ? min : x[i];
  }
  return min;
}

//
//  Masked +/=:
//
//  May be a faster way to do this with actual SIMD intrinsic masks, but just
//  putting in a placeholder for now.
//

STK_MATH_FORCE_INLINE void plus_equal_part(Double& value, const Double& increment, const int numValid) {
  assert(numValid <= ndoubles);
  for (int n=0; n<numValid; ++n) value[n] += increment[n];
}

STK_MATH_FORCE_INLINE bool are_all(const Bool& a, const int sumNum) {
  assert(sumNum <= ndoubles);
  const Double oneornone = stk::math::if_then_else_zero(a, Double(1.0));
  bool all_true = true;
  for (int i=0; i < sumNum; ++i) {
    all_true = all_true && (oneornone[i] > 0.5);
  }
  return all_true;
}

STK_MATH_FORCE_INLINE bool are_all(const Bool& a) {
#if defined(__AVX512F__) && !defined(__CUDACC__) && !defined(__HIPCC__) && !defined(USE_STK_SIMD_NONE)
  // Fast path: work directly with the __mmaskX bits.
  const unsigned k = (unsigned)a._data.get();
  const unsigned needed = ((1u << ndoubles) - 1u);
  return (k & needed) == needed;
#else
  // All lanes must be true.
  return SIMD_NAMESPACE::all_of(a._data);
#endif
}

STK_MATH_FORCE_INLINE bool are_any(const Bool& a, const int sumNum) {
  assert(sumNum <= ndoubles);
  Double oneornone = stk::math::if_then_else_zero(a, Double(1.0));
  bool any_true = false;
  for (int i=0; i < sumNum; ++i) {
    any_true = any_true || (oneornone[i] > 0.5);
  }
  return any_true;
}

STK_MATH_FORCE_INLINE bool are_any(const Bool& a) {
#if defined(__AVX512F__) && !defined(__CUDACC__) && !defined(__HIPCC__) && !defined(USE_STK_SIMD_NONE)
  const unsigned k = (unsigned)a._data.get();
  const unsigned needed = ((1u << ndoubles) - 1u);
  return (k & needed) != 0u;
#else
  return SIMD_NAMESPACE::any_of(a._data);
#endif
}

// floats

STK_MATH_FORCE_INLINE void store_part(float* x, const Float& z,const int numValid) {
  impl::specialized::store_part(x, z, numValid);
}

STK_MATH_FORCE_INLINE void store_part(float* x, const Float& z, const int offset, const int numValid) {
  impl::specialized::store_part(x, z, offset, numValid);
}

STK_MATH_FORCE_INLINE Float load_part(const float* x, const int numValid) {
  return impl::specialized::load_part(x, numValid);
}

STK_MATH_FORCE_INLINE Float load_part(const float* x, const int offset, const int numValid) {
  return impl::specialized::load_part(x, offset, numValid);
}

STK_MATH_FORCE_INLINE const float& get_data(const Float& z, int index) {
  assert(index < nfloats);
  return z[index];
}

STK_MATH_FORCE_INLINE float& get_data(Float& z, int index) {
  assert(index < nfloats);
  return z[index];
}

STK_MATH_FORCE_INLINE void set_data(Float& z, int index, const float val) {
  assert(index < nfloats);
  z[index] = val;
}

inline std::ostream& operator << (std::ostream& output, const Float& a) {
  output << "[ ";
  for ( int i=0; i<nfloats; ++i ) {
    output << a[i];
    if ( i < nfloats - 1) {
      output << ", ";
    }
  }
  output << " ]";
  return output;
}

STK_MATH_FORCE_INLINE float reduce_sum(const Float& x, const int sumNum ) {
  float ret = x[0];
  for (int i=1; i < sumNum; ++i) {
    ret += x[i];
  }
  return ret;
}

STK_MATH_FORCE_INLINE float reduce_max(const Float& x, const int sumNum ) {
  float max = x[0];
  for (int i=1; i<sumNum; ++i) {
    max = max > x[i] ? max : x[i];
  }
  return max;
}

STK_MATH_FORCE_INLINE float reduce_min(const Float& x, const int sumNum ) {
  float min = x[0];
  for (int i=1; i<sumNum; ++i) {
    min = min < x[i] ? min : x[i];
  }
  return min;
}

STK_MATH_FORCE_INLINE bool are_all(const Boolf& a, const int sumNum=nfloats) {
  assert(sumNum <= nfloats);
  const Float oneornone = stk::math::if_then_else_zero(a, Float(1.0f));
  bool all_true = true;
  for (int i=0; i < sumNum; ++i) {
    all_true = all_true && (oneornone[i] > 0.5);
  }
  return all_true;
}

STK_MATH_FORCE_INLINE bool are_any(const Boolf& a, const int sumNum=nfloats) {
  assert(sumNum <= nfloats);
  Float oneornone = stk::math::if_then_else_zero(a, Float(1.0f));
  bool any_true = false;
  for (int i=0; i < sumNum; ++i) {
    any_true = any_true || (oneornone[i] > 0.5);
  }
  return any_true;
}

STK_MATH_FORCE_INLINE int count_true(const Bool& a, const int sumNum=ndoubles) {
  assert(sumNum <= ndoubles);
  Double oneornone = stk::math::if_then_else_zero(a, Double(1.0));
  return static_cast<int>(reduce_sum(oneornone, sumNum));
}

STK_MATH_FORCE_INLINE int count_true(const Boolf& a, const int sumNum=nfloats) {
  assert(sumNum <= nfloats);
  Float oneornone = stk::math::if_then_else_zero(a, Float(1.0f));
  return static_cast<int>(reduce_sum(oneornone, sumNum));
}

STK_MATH_FORCE_INLINE const double& get_data(const double& z, [[maybe_unused]] int index) {
  assert(index==0);
  return z;
}

STK_MATH_FORCE_INLINE double& get_data(double& z, [[maybe_unused]] int index) {
  assert(index==0);
  return z;
}

STK_MATH_FORCE_INLINE void set_data(double& z, [[maybe_unused]] int index, const double val) {
  assert(index==0);
  z = val;
}

STK_MATH_FORCE_INLINE const float& get_data(const float& z, [[maybe_unused]] int index) {
  assert(index==0);
  return z;
}

STK_MATH_FORCE_INLINE float& get_data(float& z, [[maybe_unused]] int index) {
  assert(index==0);
  return z;
}

STK_MATH_FORCE_INLINE void set_data(float& z, [[maybe_unused]] int index, const float val) {
  assert(index==0);
  z = val;
}

// horizonal operations defined for scalar types for portability

STK_MATH_FORCE_INLINE double reduce_sum(const double& x) {
  return x;
}

STK_MATH_FORCE_INLINE double reduce_sum(const double& x, const int) {
  return x;
}


STK_MATH_FORCE_INLINE void plus_equal_part(double& value, const double& increment, const int /*numValid*/) {
  value += increment;
}

STK_MATH_FORCE_INLINE float reduce_sum(const float& x) {
  return x;
}

STK_MATH_FORCE_INLINE float reduce_sum(const float& x, const int) {
  return x;
}

STK_MATH_FORCE_INLINE int count_true(const bool& x, [[maybe_unused]] const int sumNum=1) {
  return x ? 1 : 0;
}

STK_MATH_FORCE_INLINE bool are_all(bool a, [[maybe_unused]] const int sumNum=1) {
  assert(sumNum==1);
  return a;
}

STK_MATH_FORCE_INLINE bool are_any(bool a, [[maybe_unused]] const int sumNum=1) {
  assert(sumNum==1);
  return a;
}

// Simd Casts: templated based on the primitive type

template <typename PRIMITIVE>
KOKKOS_FORCEINLINE_FUNCTION
typename Traits<PRIMITIVE>::simd_type* simd_ptr_cast(PRIMITIVE* x) {
  return reinterpret_cast<typename Traits<PRIMITIVE>::simd_type*>(x);
}

template <typename PRIMITIVE>
KOKKOS_FORCEINLINE_FUNCTION
typename Traits<PRIMITIVE>::simd_type const* simd_ptr_cast(PRIMITIVE const* x) {
  return reinterpret_cast<typename Traits<PRIMITIVE>::simd_type const*>(x);
}

template <typename PRIMITIVE>
KOKKOS_FORCEINLINE_FUNCTION
typename Traits<PRIMITIVE>::simd_type& simd_ref_cast(PRIMITIVE& x) {
  return *simd_ptr_cast(&x);
}

template <typename PRIMITIVE>
KOKKOS_FORCEINLINE_FUNCTION
const typename Traits<PRIMITIVE>::simd_type& simd_ref_cast(const PRIMITIVE& x) {
  return *simd_ptr_cast(&x);
}

// double versions
template <int size>
STK_MATH_FORCE_INLINE
void load_array(Double* const to, const double* const from) {
  impl::specialized::load_array<size>(to, from);
}

template <int size>
STK_MATH_FORCE_INLINE
void store_array(double* const to, const Double* const from) {
  impl::specialized::store_array<size>(to, from);
}

template <int size>
STK_MATH_FORCE_INLINE
void load_array(Double* const to, const double* const from, const int numValid) {
  impl::specialized::load_array<size>(to, from, numValid);
}

template <int size>
STK_MATH_FORCE_INLINE
void store_array(double* const to, const Double* const from, const int numValid) {
  impl::specialized::store_array<size>(to, from, numValid);
}

// float versions

template <int size>
STK_MATH_FORCE_INLINE
void load_array(Float* const to, const float* const from) {
  impl::specialized::load_array<size>(to, from);
}

template<typename T, int size>
STK_MATH_FORCE_INLINE
void store_array(float* const to, T(&from)[size]) {
  impl::specialized::store_array<size>(to, from);
}

template <int size>
STK_MATH_FORCE_INLINE
void store_array(float* const to, const Float* const from) {
  impl::specialized::store_array<size>(to, from);
}

template <int size>
STK_MATH_FORCE_INLINE
void load_array(Float* const to, const float* const from, const int numValid) {
  impl::specialized::load_array<size>(to, from, numValid);
}

template <int size>
STK_MATH_FORCE_INLINE
void store_array(float* const to, const Float* const from, const int numValid) {
  impl::specialized::store_array<size>(to, from, numValid);
}

} // namespace simd
} // namespace stk

#undef STK_INCLUDE_ONLY_STK_SIMD_HEADER
#undef STK_USE_AVX512_SPECIALIZATION

#endif // #ifndef STK_SIMD_SIMD_HPP

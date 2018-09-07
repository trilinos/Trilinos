// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_SIMD_FUNCTIONS_H
#define STK_SIMD_FUNCTIONS_H

#include "stk_simd/SimdConfig.hpp" // IWYU pragma: export
#include <iostream>

#include <stk_math/StkMath.hpp>

#define STK_HAVE_SIMD

#if defined ( STK_SIMD_AVX512 )
#include "avx512/Avx512.hpp"
#elif defined ( STK_SIMD_AVX )
#include "avx/Avx.hpp"
#elif defined ( STK_SIMD_SSE )
#include "sse/Sse.hpp"
#else
#include "no_simd/NoSimd.hpp"
#endif // Check SIMD version

#include "AlignedAllocator.hpp"
#include "Traits.hpp" // has to be included after Double, Bool, Float, Boolf are defined

#include <Kokkos_Macros.hpp>
#include <sys/time.h>
#include <sys/resource.h>

namespace stk {

inline double get_time_in_seconds() {
  timeval tp;
  struct timezone tz;
  gettimeofday(&tp, &tz);
  return (tp.tv_sec + tp.tv_usec/1000000.0);
}

namespace simd {

//
//  These are helper functions for processing remainder loops.  The end pieces of SIMD loops
//  that do not have a full 'stk::ndoubles' entries.  'numValid' are the number of actual
//  entries to be processed.  Assumed that numValid is at least one and is less than
//  stk::ndoubles
//

// double versions

STK_MATH_FORCE_INLINE void store_part(double* x, const Double& z, const int numValid) {
  assert(numValid <= ndoubles);
  for(int n=0; n<numValid; ++n) x[n] = z[n];
}

STK_MATH_FORCE_INLINE void store_part(double* x, const Double& z, const int offset, const int numValid) {
  assert(numValid <= ndoubles);
  for(int n=0; n<numValid; ++n) x[n*offset] = z[n];
}

STK_MATH_FORCE_INLINE Double load_part(const double* x, const int numValid) {
  assert(numValid <= ndoubles);
  Double temp(0.0);
  for(int n=0; n<numValid; ++n) temp[n] = x[n];
  return temp;
}

STK_MATH_FORCE_INLINE Double load_part(const double* x, const int offset, const int numValid) {
  assert(numValid <= ndoubles);
  Double temp(0.0);
  for(int n=0; n<numValid; ++n) temp[n] = x[n*offset];
  return temp;
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

//
//  Masked +/=:
//
//  May be a faster way to do this with actual SIMD intrinsic masks, but just
//  putting in a placeholder for now.
//

STK_MATH_FORCE_INLINE void plus_equal_part(Double& value, const Double& increment, const int numValid) {
  assert(numValid <= ndoubles);
  for(int n=0; n<numValid; ++n) value[n] += increment[n];
}

STK_MATH_FORCE_INLINE bool are_all(const Bool& a, const int sumNum=ndoubles) {
  assert(sumNum <= ndoubles);
  const Double oneornone = stk::math::if_then_else_zero(a, Double(1.0));
  bool all_true = true;
  for (int i=0; i < sumNum; ++i) {
    all_true = all_true && (oneornone[i] > 0.5);
  }
  return all_true;
}

STK_MATH_FORCE_INLINE bool are_any(const Bool& a, const int sumNum=ndoubles) {
  assert(sumNum <= ndoubles);
  Double oneornone = stk::math::if_then_else_zero(a, Double(1.0));
  bool any_true = false;
  for (int i=0; i < sumNum; ++i) {
    any_true = any_true || (oneornone[i] > 0.5);
  }
  return any_true;
}

// floats

STK_MATH_FORCE_INLINE void store_part(float* x, const Float& z,const int numValid) {
  assert(numValid <= nfloats);
  for(int n=0; n<numValid; ++n) x[n] = z[n];
}

STK_MATH_FORCE_INLINE void store_part(float* x, const Float& z, const int offset, const int numValid) {
  assert(numValid <= nfloats);
  for(int n=0; n<numValid; ++n) x[n*offset] = z[n];
}

STK_MATH_FORCE_INLINE Float load_part(const float* x, const int numValid) {
  assert(numValid <= nfloats);
  Float temp;
  for(int n=0; n<numValid; ++n) temp[n] = x[n];
  return temp;
}

STK_MATH_FORCE_INLINE Float load_part(const float* x, const int offset, const int numValid) {
  assert(numValid <= nfloats);
  Float temp;
  for(int n=0; n<numValid; ++n) temp[n] = x[n*offset];
  return temp;
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

STK_MATH_FORCE_INLINE const double& get_data(const double& z, int index) {
  assert(index==0);
  return z;
}

STK_MATH_FORCE_INLINE double& get_data(double& z, int index) {
  assert(index==0);
  return z;
}

STK_MATH_FORCE_INLINE void set_data(double& z, int index, const double val) {
  assert(index==0);
  z = val;
}

STK_MATH_FORCE_INLINE const float& get_data(const float& z, int index) {
  assert(index==0);
  return z;
}

STK_MATH_FORCE_INLINE float& get_data(float& z, int index) {
  assert(index==0);
  return z;
}

STK_MATH_FORCE_INLINE void set_data(float& z, int index, const float val) {
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


STK_MATH_FORCE_INLINE void plus_equal_part(double& value, const double& increment, const int numValid) {
  value += increment;
}

STK_MATH_FORCE_INLINE float reduce_sum(const float& x) {
  return x;
}

STK_MATH_FORCE_INLINE float reduce_sum(const float& x, const int) {
  return x;
}

STK_MATH_FORCE_INLINE int count_true(const bool& x, const int sumNum=1) {
  return x ? 1 : 0;
}

STK_MATH_FORCE_INLINE bool are_all(bool a, const int sumNum=1) {
  assert(sumNum==1);
  return a;
}

STK_MATH_FORCE_INLINE bool are_any(bool a, const int sumNum=1) {
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
  for (int i=0; i < size; ++i) {
    to[i] = load(from+i,size);
  }
}

template <int size>
STK_MATH_FORCE_INLINE
void store_array(double* const to, const Double* const from) {
  for (int i=0; i < size; ++i) {
    store(to+i,from[i],size);
  }
}

template <int size>
STK_MATH_FORCE_INLINE
void load_array(Double* const to, const double* const from, const int numValid) {
  for (int i=0; i < size; ++i) {
    to[i] = load_part(from+i,size,numValid);
  }
}

template <int size>
STK_MATH_FORCE_INLINE
void store_array(double* const to, const Double* const from, const int numValid) {
  for (int i=0; i < size; ++i) {
    store_part(to+i,from[i],size,numValid);
  }
}

// float versions

template <int size>
STK_MATH_FORCE_INLINE
void load_array(Float* const to, const float* const from) {
  for (int i=0; i < size; ++i) {
    to[i] = load(from+i,size);
  }
}

template<typename T, int size>
STK_MATH_FORCE_INLINE
void store_array(float* const to, T(&from)[size]) {
  for (int i=0; i < size; ++i) {
    store(to+i,from[i],size);
  }
}

template <int size>
STK_MATH_FORCE_INLINE
void store_array(float* const to, const Float* const from) {
  for (int i=0; i < size; ++i) {
    store(to+i,from[i],size);
  }
}

template <int size>
STK_MATH_FORCE_INLINE
void load_array(Float* const to, const float* const from, const int numValid) {
  for (int i=0; i < size; ++i) {
    to[i] = load_part(from+i,size,numValid);
  }
}

template <int size>
STK_MATH_FORCE_INLINE
void store_array(float* const to, const Float* const from, const int numValid) {
  for (int i=0; i < size; ++i) {
    store_part(to+i,from[i],size,numValid);
  }
}

} // namespace simd
} // namespace stk

#endif // #ifndef SIMD_H__

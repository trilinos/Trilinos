//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#ifndef __KOKKOSBATCHED_VECTOR_SIMD_HPP__
#define __KOKKOSBATCHED_VECTOR_SIMD_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <Kokkos_Complex.hpp>
#include <KokkosBatched_Vector.hpp>
#include "KokkosKernels_Macros.hpp"

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
#undef __KOKKOSBATCHED_ENABLE_AVX__
#else
// compiler bug with AVX in some architectures
#define __KOKKOSBATCHED_ENABLE_AVX__
#endif

namespace KokkosBatched {

template <typename T, int l>
class Vector<SIMD<T>, l> {
 public:
  using type       = Vector<SIMD<T>, l>;
  using value_type = T;
  using mag_type   = typename Kokkos::ArithTraits<T>::mag_type;

  enum : int { vector_length = l };

  typedef value_type data_type[vector_length];

  KOKKOS_INLINE_FUNCTION
  static const char *label() { return "SIMD"; }

  template <typename, int>
  friend class Vector;

 private:
  mutable data_type _data;

 public:
  KOKKOS_INLINE_FUNCTION Vector() {
    // NOTE Not meant to be instantiated for CUDA
    KOKKOSKERNELS_FORCE_SIMD
    for (int i = 0; i < vector_length; ++i) _data[i] = 0;
  }
  template <typename ArgValueType>
  KOKKOS_INLINE_FUNCTION Vector(const ArgValueType &val) {
    KOKKOSKERNELS_FORCE_SIMD
    for (int i = 0; i < vector_length; ++i) _data[i] = val;
  }
  template <typename ArgValueType>
  KOKKOS_INLINE_FUNCTION Vector(const Vector<SIMD<ArgValueType>, vector_length> &b) {
    KOKKOSKERNELS_FORCE_SIMD
    for (int i = 0; i < vector_length; ++i) _data[i] = b[i];
  }

  KOKKOS_INLINE_FUNCTION
  type &loadAligned(const value_type *p) {
    KOKKOSKERNELS_FORCE_SIMD
    for (int i = 0; i < vector_length; ++i) _data[i] = p[i];
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  type &loadUnaligned(const value_type *p) { return loadAligned(p); }

  KOKKOS_INLINE_FUNCTION
  void storeAligned(value_type *p) const {
    KOKKOSKERNELS_FORCE_SIMD
    for (int i = 0; i < vector_length; ++i) p[i] = _data[i];
  }

  KOKKOS_INLINE_FUNCTION
  void storeUnaligned(value_type *p) const { storeAligned(p); }

  KOKKOS_INLINE_FUNCTION
  value_type &operator[](const int &i) const { return _data[i]; }
};
}  // namespace KokkosBatched

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
namespace KokkosBatched {

template <>
class Vector<SIMD<float>, 2> {
 public:
  using type       = Vector<SIMD<float>, 2>;
  using value_type = float;
  using mag_type   = float;

  enum : int { vector_length = 2 };
  typedef float2 data_type;

  KOKKOS_INLINE_FUNCTION
  static const char *label() { return "GpuFloat2"; }

  template <typename, int>
  friend class Vector;

 private:
  mutable data_type _data;

 public:
  KOKKOS_INLINE_FUNCTION Vector() {
    _data.x = 0;
    _data.y = 0;
  }
  KOKKOS_INLINE_FUNCTION Vector(const value_type &val) {
    _data.x = val;
    _data.y = val;
  }
  KOKKOS_INLINE_FUNCTION Vector(const type &b) {
    _data.x = b._data.x;
    _data.y = b._data.y;
  }
  KOKKOS_INLINE_FUNCTION Vector(const float2 &val) {
    _data.x = val.x;
    _data.y = val.y;
  }

  template <typename ArgValueType>
  KOKKOS_INLINE_FUNCTION Vector(const ArgValueType &val) {
    _data.x = val;
    _data.y = val;
  }

  template <typename ArgValueType>
  KOKKOS_INLINE_FUNCTION Vector(const Vector<SIMD<ArgValueType>, vector_length> &b) {
    _data.x = b[0];
    _data.y = b[1];
  }

  KOKKOS_INLINE_FUNCTION
  type &operator=(const float2 &val) {
    _data.x = val.x;
    _data.y = val.y;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  float2 float2() const { return _data; }

  KOKKOS_INLINE_FUNCTION
  type &loadAligned(const value_type *p) {
    _data.x = *(p);
    _data.y = *(p + 1);
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  type &loadUnaligned(const value_type *p) {
    _data.x = *(p);
    _data.y = *(p + 1);
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void storeAligned(value_type *p) const {
    *(p)     = _data.x;
    *(p + 1) = _data.y;
  }

  KOKKOS_INLINE_FUNCTION
  void storeUnaligned(value_type *p) const {
    *(p)     = _data.x;
    *(p + 1) = _data.y;
  }

  KOKKOS_INLINE_FUNCTION
  value_type &operator[](const int &i) const { return reinterpret_cast<value_type *>(&_data)[i]; }
};

template <>
class Vector<SIMD<double>, 2> {
 public:
  using type       = Vector<SIMD<double>, 2>;
  using value_type = double;
  using mag_type   = double;

  enum : int { vector_length = 2 };
  typedef double2 data_type;

  KOKKOS_INLINE_FUNCTION
  static const char *label() { return "GpuDouble2"; }

  template <typename, int>
  friend class Vector;

 private:
  mutable data_type _data;

 public:
  KOKKOS_INLINE_FUNCTION Vector() {
    _data.x = 0;
    _data.y = 0;
  }
  KOKKOS_INLINE_FUNCTION Vector(const value_type &val) {
    _data.x = val;
    _data.y = val;
  }
  KOKKOS_INLINE_FUNCTION Vector(const type &b) {
    _data.x = b._data.x;
    _data.y = b._data.y;
  }
  KOKKOS_INLINE_FUNCTION Vector(const double2 &val) {
    _data.x = val.x;
    _data.y = val.y;
  }

  template <typename ArgValueType>
  KOKKOS_INLINE_FUNCTION Vector(const ArgValueType &val) {
    _data.x = val;
    _data.y = val;
  }

  template <typename ArgValueType>
  KOKKOS_INLINE_FUNCTION Vector(const Vector<SIMD<ArgValueType>, vector_length> &b) {
    _data.x = b[0];
    _data.y = b[1];
  }

  KOKKOS_INLINE_FUNCTION
  type &operator=(const double2 &val) {
    _data.x = val.x;
    _data.y = val.y;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  double2 double2() const { return _data; }

  KOKKOS_INLINE_FUNCTION
  type &loadAligned(const value_type *p) {
    _data.x = *(p);
    _data.y = *(p + 1);
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  type &loadUnaligned(const value_type *p) {
    _data.x = *(p);
    _data.y = *(p + 1);
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void storeAligned(value_type *p) const {
    *(p)     = _data.x;
    *(p + 1) = _data.y;
  }

  KOKKOS_INLINE_FUNCTION
  void storeUnaligned(value_type *p) const {
    *(p)     = _data.x;
    *(p + 1) = _data.y;
  }

  KOKKOS_INLINE_FUNCTION
  value_type &operator[](const int &i) const { return reinterpret_cast<value_type *>(&_data)[i]; }
};

template <>
class Vector<SIMD<float>, 4> {
 public:
  using type       = Vector<SIMD<float>, 4>;
  using value_type = float;
  using mag_type   = float;

  enum : int { vector_length = 4 };
  typedef float4 data_type;

  KOKKOS_INLINE_FUNCTION
  static const char *label() { return "GpuFloat4"; }

  template <typename, int>
  friend class Vector;

 private:
  mutable data_type _data;

 public:
  KOKKOS_INLINE_FUNCTION Vector() {
    _data.x = 0;
    _data.y = 0;
    _data.z = 0;
    _data.w = 0;
  }
  KOKKOS_INLINE_FUNCTION Vector(const value_type &val) {
    _data.x = val;
    _data.y = val;
    _data.z = val;
    _data.w = val;
  }
  KOKKOS_INLINE_FUNCTION Vector(const type &b) {
    _data.x = b._data.x;
    _data.y = b._data.y;
    _data.z = b._data.z;
    _data.w = b._data.w;
  }
  KOKKOS_INLINE_FUNCTION Vector(const float4 &val) {
    _data.x = val.x;
    _data.y = val.y;
    _data.z = val.z;
    _data.w = val.w;
  }

  template <typename ArgValueType>
  KOKKOS_INLINE_FUNCTION Vector(const ArgValueType &val) {
    _data.x = val;
    _data.y = val;
    _data.z = val;
    _data.w = val;
  }

  template <typename ArgValueType>
  KOKKOS_INLINE_FUNCTION Vector(const Vector<SIMD<ArgValueType>, vector_length> &b) {
    _data.x = b[0];
    _data.y = b[1];
    _data.z = b[2];
    _data.w = b[3];
  }

  KOKKOS_INLINE_FUNCTION
  type &operator=(const float4 &val) {
    _data.x = val.x;
    _data.y = val.y;
    _data.z = val.z;
    _data.w = val.w;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  float4 float4() const { return _data; }

  KOKKOS_INLINE_FUNCTION
  type &loadAligned(const value_type *p) {
    _data.x = *(p);
    _data.y = *(p + 1);
    _data.z = *(p + 2);
    _data.w = *(p + 3);
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  type &loadUnaligned(const value_type *p) {
    _data.x = *(p);
    _data.y = *(p + 1);
    _data.z = *(p + 2);
    _data.w = *(p + 3);
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void storeAligned(value_type *p) const {
    *(p)     = _data.x;
    *(p + 1) = _data.y;
    *(p + 2) = _data.z;
    *(p + 3) = _data.w;
  }

  KOKKOS_INLINE_FUNCTION
  void storeUnaligned(value_type *p) const {
    *(p)     = _data.x;
    *(p + 1) = _data.y;
    *(p + 2) = _data.z;
    *(p + 3) = _data.w;
  }

  KOKKOS_INLINE_FUNCTION
  value_type &operator[](const int &i) const { return reinterpret_cast<value_type *>(&_data)[i]; }
};

template <>
class Vector<SIMD<double>, 4> {
 public:
  using type       = Vector<SIMD<double>, 4>;
  using value_type = double;
  using mag_type   = double;

  enum : int { vector_length = 4 };
  typedef double4 data_type;

  KOKKOS_INLINE_FUNCTION
  static const char *label() { return "GpuDouble4"; }

  template <typename, int>
  friend class Vector;

 private:
  mutable data_type _data;

 public:
  KOKKOS_INLINE_FUNCTION Vector() {
    _data.x = 0;
    _data.y = 0;
    _data.z = 0;
    _data.w = 0;
  }
  KOKKOS_INLINE_FUNCTION Vector(const value_type &val) {
    _data.x = val;
    _data.y = val;
    _data.z = val;
    _data.w = val;
  }
  KOKKOS_INLINE_FUNCTION Vector(const type &b) {
    _data.x = b._data.x;
    _data.y = b._data.y;
    _data.z = b._data.z;
    _data.w = b._data.w;
  }
  KOKKOS_INLINE_FUNCTION Vector(const double4 &val) {
    _data.x = val.x;
    _data.y = val.y;
    _data.z = val.z;
    _data.w = val.w;
  }

  template <typename ArgValueType>
  KOKKOS_INLINE_FUNCTION Vector(const ArgValueType &val) {
    _data.x = val;
    _data.y = val;
    _data.z = val;
    _data.w = val;
  }

  template <typename ArgValueType>
  KOKKOS_INLINE_FUNCTION Vector(const Vector<SIMD<ArgValueType>, vector_length> &b) {
    _data.x = b[0];
    _data.y = b[1];
    _data.z = b[2];
    _data.w = b[3];
  }

  KOKKOS_INLINE_FUNCTION
  type &operator=(const double4 &val) {
    _data.x = val.x;
    _data.y = val.y;
    _data.z = val.z;
    _data.w = val.w;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  double4 double4() const { return _data; }

  KOKKOS_INLINE_FUNCTION
  type &loadAligned(const value_type *p) {
    _data.x = *(p);
    _data.y = *(p + 1);
    _data.z = *(p + 2);
    _data.w = *(p + 3);
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  type &loadUnaligned(const value_type *p) {
    _data.x = *(p);
    _data.y = *(p + 1);
    _data.z = *(p + 2);
    _data.w = *(p + 3);
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void storeAligned(value_type *p) const {
    *(p)     = _data.x;
    *(p + 1) = _data.y;
    *(p + 2) = _data.z;
    *(p + 3) = _data.w;
  }

  KOKKOS_INLINE_FUNCTION
  void storeUnaligned(value_type *p) const {
    *(p)     = _data.x;
    *(p + 1) = _data.y;
    *(p + 2) = _data.z;
    *(p + 3) = _data.w;
  }

  KOKKOS_INLINE_FUNCTION
  value_type &operator[](const int &i) const { return reinterpret_cast<value_type *>(&_data)[i]; }
};

}  // namespace KokkosBatched
#endif

#if defined(__KOKKOSBATCHED_ENABLE_AVX__)
#if defined(__AVX__) || defined(__AVX2__)

#if CUDA_VERSION < 12022
#undef _Float16
#endif

#include <immintrin.h>

namespace KokkosBatched {

template <>
class Vector<SIMD<double>, 4> {
 public:
  using type       = Vector<SIMD<double>, 4>;
  using value_type = double;
  using mag_type   = double;

  enum : int { vector_length = 4 };
  typedef __m256d data_type __attribute__((aligned(32)));

  inline static const char *label() { return "AVX256"; }

  template <typename, int>
  friend class Vector;

 private:
  mutable data_type _data;

 public:
  inline Vector() { _data = _mm256_setzero_pd(); }
  inline Vector(const value_type &val) { _data = _mm256_set1_pd(val); }
  inline Vector(const type &b) { _data = b._data; }
  inline Vector(const __m256d &val) { _data = val; }

  template <typename ArgValueType>
  inline Vector(const ArgValueType &val) {
    auto d = reinterpret_cast<value_type *>(&_data);
    KOKKOSKERNELS_FORCE_SIMD
    for (int i = 0; i < vector_length; ++i) d[i] = val;
  }

  template <typename ArgValueType>
  inline Vector(const Vector<SIMD<ArgValueType>, vector_length> &b) {
    auto dd = reinterpret_cast<value_type *>(&_data);
    auto bb = reinterpret_cast<ArgValueType *>(&b._data);
    KOKKOSKERNELS_FORCE_SIMD
    for (int i = 0; i < vector_length; ++i) dd[i] = bb[i];
  }

  inline type &operator=(const __m256d &val) {
    _data = val;
    return *this;
  }

  inline operator __m256d() const { return _data; }

  inline type &loadAligned(const value_type *p) {
    _data = _mm256_load_pd(p);
    return *this;
  }

  inline type &loadUnaligned(const value_type *p) {
    _data = _mm256_loadu_pd(p);
    return *this;
  }

  inline void storeAligned(value_type *p) const { _mm256_store_pd(p, _data); }

  inline void storeUnaligned(value_type *p) const { _mm256_storeu_pd(p, _data); }

  inline value_type &operator[](const int &i) const { return reinterpret_cast<value_type *>(&_data)[i]; }
};

template <>
class Vector<SIMD<Kokkos::complex<double> >, 2> {
 public:
  using type       = Vector<SIMD<Kokkos::complex<double> >, 2>;
  using value_type = Kokkos::complex<double>;
  using mag_type   = double;

  static const int vector_length = 2;
  typedef __m256d data_type __attribute__((aligned(32)));

  inline static const char *label() { return "AVX256"; }

  template <typename, int>
  friend class Vector;

 private:
  mutable data_type _data;

 public:
  inline Vector() { _data = _mm256_setzero_pd(); }
  inline Vector(const value_type &val) {
    _data = _mm256_broadcast_pd((const __m128d *)&val);
    KOKKOSKERNELS_GNU_COMPILER_FENCE
  }
  inline Vector(const mag_type &val) {
    const value_type a(val);
    _data = _mm256_broadcast_pd((__m128d const *)&a);
    KOKKOSKERNELS_GNU_COMPILER_FENCE
  }
  inline Vector(const type &b) { _data = b._data; }
  inline Vector(const __m256d &val) { _data = val; }

  //       template<typename ArgValueType>
  //       inline Vector(const ArgValueType val) {
  // #if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
  // #pragma ivdep
  // #endif
  // #if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
  // #pragma vector always
  // #endif
  //         for (int i=0;i<vector_length;++i)
  //           _data.d[i] = value_type(val);
  //       }
  template <typename ArgValueType>
  inline Vector(const Vector<SIMD<ArgValueType>, vector_length> &b) {
    auto dd = reinterpret_cast<value_type *>(&_data);
    auto bb = reinterpret_cast<ArgValueType *>(&b._data);
    KOKKOSKERNELS_FORCE_SIMD
    for (int i = 0; i < vector_length; ++i) dd[i] = bb[i];
  }

  inline type &operator=(const __m256d &val) {
    _data = val;
    return *this;
  }

  inline operator __m256d() const { return _data; }

  inline type &loadAligned(const value_type *p) {
    _data = _mm256_load_pd((mag_type *)p);
    return *this;
  }

  inline type &loadUnaligned(const value_type *p) {
    _data = _mm256_loadu_pd((mag_type *)p);
    return *this;
  }

  inline void storeAligned(value_type *p) const { _mm256_store_pd((mag_type *)p, _data); }

  inline void storeUnaligned(value_type *p) const { _mm256_storeu_pd((mag_type *)p, _data); }

  inline value_type &operator[](const int &i) const { return reinterpret_cast<value_type *>(&_data)[i]; }
};
}  // namespace KokkosBatched
#endif /* #if defined(__AVX__) || defined(__AVX2__) */

#if defined(__AVX512F__)
#if CUDA_VERSION < 12022
#undef _Float16
#endif
#include <immintrin.h>

namespace KokkosBatched {

template <>
class Vector<SIMD<double>, 8> {
 public:
  using type       = Vector<SIMD<double>, 8>;
  using value_type = double;
  using mag_type   = double;

  enum : int { vector_length = 8 };
  typedef __m512d data_type __attribute__((aligned(64)));

  inline static const char *label() { return "AVX512"; }

  template <typename, int>
  friend class Vector;

 private:
  mutable data_type _data;

 public:
  inline Vector() { _data = _mm512_setzero_pd(); }
  inline Vector(const value_type &val) { _data = _mm512_set1_pd(val); }
  inline Vector(const type &b) { _data = b._data; }
  inline Vector(const __m512d &val) { _data = val; }

  template <typename ArgValueType>
  inline Vector(const ArgValueType &val) {
    auto d = reinterpret_cast<value_type *>(&_data);
    KOKKOSKERNELS_FORCE_SIMD
    for (int i = 0; i < vector_length; ++i) d[i] = val;
  }
  template <typename ArgValueType>
  inline Vector(const Vector<SIMD<ArgValueType>, vector_length> &b) {
    auto dd = reinterpret_cast<value_type *>(&_data);
    auto bb = reinterpret_cast<ArgValueType *>(&b._data);
    KOKKOSKERNELS_FORCE_SIMD
    for (int i = 0; i < vector_length; ++i) dd[i] = bb[i];
  }

  inline type &operator=(const __m512d &val) {
    _data = val;
    return *this;
  }

  inline operator __m512d() const { return _data; }

  inline type &loadAligned(const value_type *p) {
    _data = _mm512_load_pd(p);
    return *this;
  }

  inline type &loadUnaligned(const value_type *p) {
    _data = _mm512_loadu_pd(p);
    return *this;
  }

  inline void storeAligned(value_type *p) const { _mm512_store_pd(p, _data); }

  inline void storeUnaligned(value_type *p) const { _mm512_storeu_pd(p, _data); }

  inline value_type &operator[](const int &i) const { return reinterpret_cast<value_type *>(&_data)[i]; }
};

template <>
class Vector<SIMD<Kokkos::complex<double> >, 4> {
 public:
  using type       = Vector<SIMD<Kokkos::complex<double> >, 4>;
  using value_type = Kokkos::complex<double>;
  using mag_type   = double;

  enum : int { vector_length = 4 };
  typedef __m512d data_type __attribute__((aligned(64)));

  inline static const char *label() { return "AVX512"; }

  template <typename, int>
  friend class Vector;

 private:
  mutable data_type _data;

 public:
  inline Vector() { _data = _mm512_setzero_pd(); }
  inline Vector(const value_type &val) {
    _data = _mm512_mask_broadcast_f64x4(_mm512_set1_pd(val.imag()), 0x55, _mm256_set1_pd(val.real()));
    KOKKOSKERNELS_GNU_COMPILER_FENCE
  }
  inline Vector(const mag_type &val) {
    _data = _mm512_mask_broadcast_f64x4(_mm512_setzero_pd(), 0x55, _mm256_set1_pd(val));
    KOKKOSKERNELS_GNU_COMPILER_FENCE
  }
  inline Vector(const type &b) { _data = b._data; }
  inline Vector(const __m512d &val) { _data = val; }

  template <typename ArgValueType>
  inline Vector(const ArgValueType &val) {
    auto d = reinterpret_cast<value_type *>(&_data);
    KOKKOSKERNELS_FORCE_SIMD
    for (int i = 0; i < vector_length; ++i) d[i] = val;
  }
  template <typename ArgValueType>
  inline Vector(const Vector<SIMD<ArgValueType>, vector_length> &b) {
    auto dd = reinterpret_cast<value_type *>(&_data);
    auto bb = reinterpret_cast<value_type *>(&b._data);
    KOKKOSKERNELS_FORCE_SIMD
    for (int i = 0; i < vector_length; ++i) dd[i] = bb[i];
  }

  inline type &operator=(const __m512d &val) {
    _data = val;
    return *this;
  }

  inline operator __m512d() const { return _data; }

  inline type &loadAligned(const value_type *p) {
    _data = _mm512_load_pd((mag_type *)p);
    return *this;
  }

  inline type &loadUnaligned(const value_type *p) {
    _data = _mm512_loadu_pd((mag_type *)p);
    return *this;
  }

  inline void storeAligned(value_type *p) const { _mm512_store_pd((mag_type *)p, _data); }

  inline void storeUnaligned(value_type *p) const { _mm512_storeu_pd((mag_type *)p, _data); }

  inline value_type &operator[](const int &i) const { return reinterpret_cast<value_type *>(&_data)[i]; }
};
}  // namespace KokkosBatched

#endif /* #if defined(__AVX512F__) */
#endif /* #if defined(__KOKKOSBATCHED_ENABLE_AVX__) */

#include "KokkosBatched_Vector_SIMD_Arith.hpp"
#include "KokkosBatched_Vector_SIMD_Logical.hpp"
#include "KokkosBatched_Vector_SIMD_Relation.hpp"
#include "KokkosBatched_Vector_SIMD_Math.hpp"
#include "KokkosBatched_Vector_SIMD_Misc.hpp"
#include "KokkosBatched_Vector_SIMD_View.hpp"

#endif

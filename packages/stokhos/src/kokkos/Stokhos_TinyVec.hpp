// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_TINY_VEC_HPP
#define STOKHOS_TINY_VEC_HPP

#include "Stokhos_ConfigDefs.h"
#if defined(HAVE_STOKHOS_INTRINSICS) && !defined( __CUDACC__ )

extern "C" {
#include <immintrin.h>
}

#endif

#include "Kokkos_Macros.hpp"

namespace Stokhos {

#if defined(__INTEL_COMPILER) && ! defined( __CUDA_ARCH__)

template <typename ValueType, int N, bool UseIntrinsics, bool Mask = false >
class TinyVec {
public:

  static const int Num = N;

  KOKKOS_INLINE_FUNCTION
  TinyVec() {}

  KOKKOS_INLINE_FUNCTION
  TinyVec(const ValueType a[]) {
    load(a);
  }

  template <typename OrdinalType>
  KOKKOS_INLINE_FUNCTION
  TinyVec(const ValueType a[], const OrdinalType idx[]) {
    gather(a,idx);
  }

  KOKKOS_INLINE_FUNCTION
  TinyVec(const ValueType a) {
    load(a);
  }

  KOKKOS_INLINE_FUNCTION
  TinyVec(const TinyVec& tv) {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<Num; ++i)
      v[i] = tv.v[i];
  }

  KOKKOS_INLINE_FUNCTION
  TinyVec& operator=(const TinyVec& tv) {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<Num; ++i)
      v[i] = tv.v[i];
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void load(const ValueType a[]) {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<Num; ++i)
      v[i] = a[i];
  }

  KOKKOS_INLINE_FUNCTION
  void load(const ValueType a) {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<Num; ++i)
      v[i] = a;
  }

  KOKKOS_INLINE_FUNCTION
  void aligned_load(const ValueType a[]) {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<Num; ++i)
      v[i] = a[i];
  }

  template <typename OrdinalType>
  KOKKOS_INLINE_FUNCTION
  void gather(const ValueType a[], const OrdinalType idx[]) {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<Num; ++i)
      v[i] = a[idx[i]];
  }

  KOKKOS_INLINE_FUNCTION
  void scatter(ValueType a[]) const {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<Num; ++i)
      a[i] = v[i];
  }

  KOKKOS_INLINE_FUNCTION
  void aligned_scatter(ValueType a[]) const {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<Num; ++i)
      a[i] = v[i];
  }

  KOKKOS_INLINE_FUNCTION
  void zero() {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<Num; ++i)
      v[i] = ValueType(0.0);
  }

  KOKKOS_INLINE_FUNCTION
  void plus_equal(const TinyVec& t) {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<Num; ++i)
      v[i] += t.v[i];
  }

  KOKKOS_INLINE_FUNCTION
  void times_equal(const TinyVec& t) {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<Num; ++i)
      v[i] *= t.v[i];
  }

  // *this = *this + t1 * t2
  KOKKOS_INLINE_FUNCTION
  void multiply_add(const TinyVec& t1, const TinyVec& t2) {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<Num; ++i)
      v[i] += t1.v[i]*t2.v[i];
  }

  KOKKOS_INLINE_FUNCTION
  ValueType sum() const {
    ValueType s(0.0);
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<Num; ++i)
      s += v[i];
    return s;
  }

private:
   ValueType v[Num] __attribute__((aligned(64)));
};

template <typename ValueType, int N, bool UseIntrinsics >
class TinyVec<ValueType,N,UseIntrinsics,true> {
public:

  static const int Num = N;

  KOKKOS_INLINE_FUNCTION
  TinyVec(int size) { sz = size; }

  KOKKOS_INLINE_FUNCTION
  TinyVec(const ValueType a[], int size) {
    sz = size;
    load(a);
  }

  template <typename OrdinalType>
  KOKKOS_INLINE_FUNCTION
  TinyVec(const ValueType a[], const OrdinalType idx[], int size) {
    sz = size;
    gather(a,idx);
  }

  KOKKOS_INLINE_FUNCTION
  TinyVec(const ValueType a, int size) {
    sz = size;
    load(a);
  }

  KOKKOS_INLINE_FUNCTION
  TinyVec(const TinyVec& tv) {
    sz = tv.sz;
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<sz; ++i)
      v[i] = tv.v[i];
  }

  KOKKOS_INLINE_FUNCTION
  TinyVec& operator=(const TinyVec& tv) {
    sz = tv.sz;
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<sz; ++i)
      v[i] = tv.v[i];
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void load(const ValueType a[]) {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<sz; ++i)
      v[i] = a[i];
  }

  KOKKOS_INLINE_FUNCTION
  void load(const ValueType a) {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<sz; ++i)
      v[i] = a;
  }

  KOKKOS_INLINE_FUNCTION
  void aligned_load(const ValueType a[]) {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<sz; ++i)
      v[i] = a[i];
  }

  template <typename OrdinalType>
  KOKKOS_INLINE_FUNCTION
  void gather(const ValueType a[], const OrdinalType idx[]) {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<sz; ++i)
      v[i] = a[idx[i]];
  }

  KOKKOS_INLINE_FUNCTION
  void scatter(ValueType a[]) const {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<sz; ++i)
      a[i] = v[i];
  }

  KOKKOS_INLINE_FUNCTION
  void aligned_scatter(ValueType a[]) const {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<sz; ++i)
      a[i] = v[i];
  }

  KOKKOS_INLINE_FUNCTION
  void zero() {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<sz; ++i)
      v[i] = ValueType(0.0);
  }

  KOKKOS_INLINE_FUNCTION
  void plus_equal(const TinyVec& t) {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<sz; ++i)
      v[i] += t.v[i];
  }

  KOKKOS_INLINE_FUNCTION
  void times_equal(const TinyVec& t) {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<sz; ++i)
      v[i] *= t.v[i];
  }

  // *this = *this + t1 * t2
  KOKKOS_INLINE_FUNCTION
  void multiply_add(const TinyVec& t1, const TinyVec& t2) {
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<sz; ++i)
      v[i] += t1.v[i]*t2.v[i];
  }

  KOKKOS_INLINE_FUNCTION
  ValueType sum() const {
    ValueType s(0.0);
#pragma ivdep
#pragma vector aligned
    for (int i=0; i<sz; ++i)
      s += v[i];
    return s;
  }

private:
  ValueType v[Num] __attribute__((aligned(64)));
  int sz;
};

#else

template <typename ValueType, int N, bool UseIntrinsics, bool Mask = false >
class TinyVec {
public:

  static const int Num = N;

  KOKKOS_INLINE_FUNCTION
  TinyVec() {}

  KOKKOS_INLINE_FUNCTION
  TinyVec(const ValueType a[]) {
    load(a);
  }

  template <typename OrdinalType>
  KOKKOS_INLINE_FUNCTION
  TinyVec(const ValueType a[], const OrdinalType idx[]) {
    gather(a,idx);
  }

  KOKKOS_INLINE_FUNCTION
  TinyVec(const ValueType a) {
    load(a);
  }

  KOKKOS_INLINE_FUNCTION
  TinyVec(const TinyVec& tv) {
    for (int i=0; i<Num; ++i)
      v[i] = tv.v[i];
  }

  KOKKOS_INLINE_FUNCTION
  TinyVec& operator=(const TinyVec& tv) {
    for (int i=0; i<Num; ++i)
      v[i] = tv.v[i];
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void load(const ValueType a[]) {
    for (int i=0; i<Num; ++i)
      v[i] = a[i];
  }

  KOKKOS_INLINE_FUNCTION
  void load(const ValueType a) {
    for (int i=0; i<Num; ++i)
      v[i] = a;
  }

  KOKKOS_INLINE_FUNCTION
  void aligned_load(const ValueType a[]) {
    for (int i=0; i<Num; ++i)
      v[i] = a[i];
  }

  template <typename OrdinalType>
  KOKKOS_INLINE_FUNCTION
  void gather(const ValueType a[], const OrdinalType idx[]) {
    for (int i=0; i<Num; ++i)
      v[i] = a[idx[i]];
  }

  KOKKOS_INLINE_FUNCTION
  void scatter(ValueType a[]) const {
    for (int i=0; i<Num; ++i)
      a[i] = v[i];
  }

  KOKKOS_INLINE_FUNCTION
  void aligned_scatter(ValueType a[]) const {
    for (int i=0; i<Num; ++i)
      a[i] = v[i];
  }

  KOKKOS_INLINE_FUNCTION
  void zero() {
    for (int i=0; i<Num; ++i)
      v[i] = ValueType(0.0);
  }

  KOKKOS_INLINE_FUNCTION
  void plus_equal(const TinyVec& t) {
    for (int i=0; i<Num; ++i)
      v[i] += t.v[i];
  }

  KOKKOS_INLINE_FUNCTION
  void times_equal(const TinyVec& t) {
    for (int i=0; i<Num; ++i)
      v[i] *= t.v[i];
  }

  // *this = *this + t1 * t2
  KOKKOS_INLINE_FUNCTION
  void multiply_add(const TinyVec& t1, const TinyVec& t2) {
    for (int i=0; i<Num; ++i)
      v[i] += t1.v[i]*t2.v[i];
  }

  KOKKOS_INLINE_FUNCTION
  ValueType sum() const {
    ValueType s(0.0);
    for (int i=0; i<Num; ++i)
      s += v[i];
    return s;
  }

private:
   ValueType v[Num];
};

template <typename ValueType, int N, bool UseIntrinsics >
class TinyVec<ValueType,N,UseIntrinsics,true> {
public:

  static const int Num = N;

  KOKKOS_INLINE_FUNCTION
  TinyVec(int size) { sz = size; }

  KOKKOS_INLINE_FUNCTION
  TinyVec(const ValueType a[], int size) {
    sz = size;
    load(a);
  }

  template <typename OrdinalType>
  KOKKOS_INLINE_FUNCTION
  TinyVec(const ValueType a[], const OrdinalType idx[], int size) {
    sz = size;
    gather(a,idx);
  }

  KOKKOS_INLINE_FUNCTION
  TinyVec(const ValueType a, int size) {
    sz = size;
    load(a);
  }

  KOKKOS_INLINE_FUNCTION
  TinyVec(const TinyVec& tv) {
    sz = tv.sz;
    for (int i=0; i<sz; ++i)
      v[i] = tv.v[i];
  }

  KOKKOS_INLINE_FUNCTION
  TinyVec& operator=(const TinyVec& tv) {
    sz = tv.sz;
    for (int i=0; i<sz; ++i)
      v[i] = tv.v[i];
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void load(const ValueType a[]) {
    for (int i=0; i<sz; ++i)
      v[i] = a[i];
  }

  KOKKOS_INLINE_FUNCTION
  void load(const ValueType a) {
    for (int i=0; i<sz; ++i)
      v[i] = a;
  }

  KOKKOS_INLINE_FUNCTION
  void aligned_load(const ValueType a[]) {
    for (int i=0; i<sz; ++i)
      v[i] = a[i];
  }

  template <typename OrdinalType>
  KOKKOS_INLINE_FUNCTION
  void gather(const ValueType a[], const OrdinalType idx[]) {
    for (int i=0; i<sz; ++i)
      v[i] = a[idx[i]];
  }

  KOKKOS_INLINE_FUNCTION
  void scatter(ValueType a[]) const {
    for (int i=0; i<sz; ++i)
      a[i] = v[i];
  }

  KOKKOS_INLINE_FUNCTION
  void aligned_scatter(ValueType a[]) const {
    for (int i=0; i<sz; ++i)
      a[i] = v[i];
  }

  KOKKOS_INLINE_FUNCTION
  void zero() {
    for (int i=0; i<sz; ++i)
      v[i] = ValueType(0.0);
  }

  KOKKOS_INLINE_FUNCTION
  void plus_equal(const TinyVec& t) {
    for (int i=0; i<sz; ++i)
      v[i] += t.v[i];
  }

  KOKKOS_INLINE_FUNCTION
  void times_equal(const TinyVec& t) {
    for (int i=0; i<sz; ++i)
      v[i] *= t.v[i];
  }

  // *this = *this + t1 * t2
  KOKKOS_INLINE_FUNCTION
  void multiply_add(const TinyVec& t1, const TinyVec& t2) {
    for (int i=0; i<sz; ++i)
      v[i] += t1.v[i]*t2.v[i];
  }

  KOKKOS_INLINE_FUNCTION
  ValueType sum() const {
    ValueType s(0.0);
    for (int i=0; i<sz; ++i)
      s += v[i];
    return s;
  }

private:
  ValueType v[Num];
  int sz;
};

#endif

#if defined(HAVE_STOKHOS_INTRINSICS) && !defined( __CUDACC__ )

#ifdef __SSE2__
template <>
class TinyVec<double,2,true,false> {
public:

  typedef double ValueType;
  static const int Num = 2;

  TinyVec() {}

  TinyVec(const ValueType a[]) {
    load(a);
  }

  template <typename OrdinalType>
  TinyVec(const ValueType a[], const OrdinalType idx[]) {
    gather(a,idx);
  }

  TinyVec(const ValueType a) {
    load(a);
  }

  TinyVec(const TinyVec& tv) {
    v = tv.v;
  }

  TinyVec& operator=(const TinyVec& tv) {
    v = tv.v;
    return *this;
  }

  void load(const ValueType a[]) {
    v = _mm_set_pd(a[1], a[0]);
  }

  void load(const ValueType a) {
    v = _mm_set1_pd(a);
  }

  void aligned_load(const ValueType a[]) {
    v = _mm_load_pd(a);
  }

  template <typename OrdinalType>
  void gather(const ValueType a[], const OrdinalType idx[]) {
    v = _mm_set_pd(a[idx[1]], a[idx[0]]);
  }

  void scatter(ValueType a[]) const {
    _mm_storel_pd(&a[0], v);
    _mm_storeh_pd(&a[1], v);
  }

  void aligned_scatter(ValueType a[]) const {
    _mm_store_pd(a, v);
  }

  void zero() {
    v = _mm_setzero_pd();
  }

  void plus_equal(const TinyVec& t) {
    v = _mm_add_pd(v, t.v);
  }

  void times_equal(const TinyVec& t) {
    v = _mm_mul_pd(v, t.v);
  }

  // *this = *this + t1 * t2
  void multiply_add(const TinyVec& t1, const TinyVec& t2) {
    __m128d t = _mm_mul_pd(t1.v, t2.v);
    v = _mm_add_pd(v, t);
  }

  ValueType sum() const {
    ValueType a[Num];
    scatter(a);
    return a[0]+a[1];
  }

private:
  __m128d v;
};
#endif

#ifdef __AVX__
template <>
class TinyVec<float,8,true,false> {
public:

  typedef float ValueType;
  static const int Num = 8;

  TinyVec() {}

  TinyVec(const ValueType a[]) {
    aligned_load(a);
  }

  template <typename OrdinalType>
  TinyVec(const ValueType a[], const OrdinalType idx[]) {
    gather(a,idx);
  }

  TinyVec(const ValueType a) {
    load(a);
  }

  TinyVec(const TinyVec& tv) {
    v = tv.v;
  }

  TinyVec& operator=(const TinyVec& tv) {
    v = tv.v;
    return *this;
  }

  void load(const ValueType a[]) {
    v = _mm256_loadu_ps(a);
  }

  void load(const ValueType a) {
    v = _mm256_set1_ps(a);
  }

  void aligned_load(const ValueType a[]) {
    v = _mm256_load_ps(a);
  }

  template <typename OrdinalType>
  void gather(const ValueType a[], const OrdinalType idx[]) {
    __m128 v1 = _mm_set_ps(a[idx[3]], a[idx[2]], a[idx[1]], a[idx[0]]);
    __m128 v2 = _mm_set_ps(a[idx[7]], a[idx[6]], a[idx[5]], a[idx[4]]);
    v = _mm256_insertf128_ps(v, v1, 0);
    v = _mm256_insertf128_ps(v, v2, 1);
  }

  void scatter(ValueType a[]) const {
    _mm256_storeu_ps(a, v);
  }

  void aligned_scatter(ValueType a[]) const {
    _mm256_store_ps(a, v);
  }

  void zero() {
    v = _mm256_setzero_ps();
  }

  void plus_equal(const TinyVec& t) {
    v = _mm256_add_ps(v, t.v);
  }

  void times_equal(const TinyVec& t) {
    v = _mm256_mul_ps(v, t.v);
  }

  // *this = *this + t1 * t2
  void multiply_add(const TinyVec& t1, const TinyVec& t2) {
    __m256 t = _mm256_mul_ps(t1.v, t2.v);
    v = _mm256_add_ps(v, t);
  }

  ValueType sum() {
    __m256 s = _mm256_hadd_ps(v,v);
    __m128 sl = _mm256_extractf128_ps(s, 0);
    __m128 sh = _mm256_extractf128_ps(s, 1);
    sl = _mm_add_ps(sl,sh);
    sl = _mm_hadd_ps(sl,sl);
    ValueType res;
    _MM_EXTRACT_FLOAT(res, sl, 0);

    return res;
  }

private:
  __m256 v;
};

template <>
class TinyVec<double,4,true,false> {
public:

  typedef double ValueType;
  static const int Num = 4;

  TinyVec() {}

  TinyVec(const ValueType a[]) {
    aligned_load(a);
  }

  template <typename OrdinalType>
  TinyVec(const ValueType a[], const OrdinalType idx[]) {
    gather(a,idx);
  }

  TinyVec(const ValueType a) {
    load(a);
  }

  TinyVec(const TinyVec& tv) {
    v = tv.v;
  }

  TinyVec& operator=(const TinyVec& tv) {
    v = tv.v;
    return *this;
  }

  void load(const ValueType a[]) {
    v = _mm256_loadu_pd(a);
  }

  void load(const ValueType a) {
    v = _mm256_set1_pd(a);
  }

  void aligned_load(const ValueType a[]) {
    v = _mm256_load_pd(a);
  }

  template <typename OrdinalType>
  void gather(const ValueType a[], const OrdinalType idx[]) {
    __m128d v1 = _mm_set_pd(a[idx[1]], a[idx[0]]);
    __m128d v2 = _mm_set_pd(a[idx[3]], a[idx[2]]);
    v = _mm256_insertf128_pd(v, v1, 0);
    v = _mm256_insertf128_pd(v, v2, 1);
  }

  void scatter(ValueType a[]) const {
    _mm256_storeu_pd(a, v);
  }

  void aligned_scatter(ValueType a[]) const {
    _mm256_store_pd(a, v);
  }

  void zero() {
    v = _mm256_setzero_pd();
  }

  void plus_equal(const TinyVec& t) {
    v = _mm256_add_pd(v, t.v);
  }

  void times_equal(const TinyVec& t) {
    v = _mm256_mul_pd(v, t.v);
  }

  // *this = *this + t1 * t2
  void multiply_add(const TinyVec& t1, const TinyVec& t2) {
    __m256d t = _mm256_mul_pd(t1.v, t2.v);
    v = _mm256_add_pd(v, t);
  }

  ValueType sum() {
    // ValueType a[Num];
    // scatter(a);
    // return a[0]+a[1]+a[2]+a[3];

    // __m128d vl = _mm256_extractf128_pd(v, 0); // v[0], v[1]
    // __m128d vh = _mm256_extractf128_pd(v, 1); // v[2], v[3]
    // vh = _mm_hadd_pd(vl, vh); // v[0]+v[1], v[2]+v[3]
    // vh = _mm_hadd_pd(vh, vh); // v[0]+v[1]+v[2]+v[3], v[0]+v[1]+v[2]+v[3]
    // ValueType res;
    // _mm_storel_pd(&res, vh);
    // return res;

    __m256d s = _mm256_hadd_pd(v,v); //v[0]+v[1] v[0]+v[1] v[2]+v[3] v[2]+v[3]
    __m128d sl = _mm256_extractf128_pd(s, 0); //v[0]+v[1] v[0]+v[1]
    __m128d sh = _mm256_extractf128_pd(s, 1); //v[2]+v[3] v[2]+v[3]
    sl = _mm_add_pd(sl,sh); // v[0]+v[1]+v[2]+v[3] v[0]+v[1]+v[2]+v[3]
    ValueType res;
    _mm_storel_pd(&res, sl);
    return res;
  }

private:
  __m256d v;
};

template <>
class TinyVec<double,8,true,false> {
public:

  typedef double ValueType;
  static const int Num = 8;

  TinyVec() {}

  TinyVec(const ValueType a[]) {
    load(a);
  }

  template <typename OrdinalType>
  TinyVec(const ValueType a[], const OrdinalType idx[]) {
    gather(a,idx);
  }

  TinyVec(const ValueType a) {
    load(a);
  }

  TinyVec(const TinyVec& tv) {
    v1 = tv.v1; v2 = tv.v2;
  }

  TinyVec& operator=(const TinyVec& tv) {
    v1 = tv.v1; v2 = tv.v2;
    return *this;
  }

  void load(const ValueType a[]) {
    v1 = _mm256_loadu_pd(a);
    v2 = _mm256_loadu_pd(a+4);
  }

  void load(const ValueType a) {
    v1 = _mm256_set1_pd(a);
    v2 = _mm256_set1_pd(a);
  }

  void aligned_load(const ValueType a[]) {
    v1 = _mm256_load_pd(a);
    v2 = _mm256_load_pd(a+4);
  }

  template <typename OrdinalType>
  void gather(const ValueType a[], const OrdinalType idx[]) {
    __m128d t1 = _mm_set_pd(a[idx[1]], a[idx[0]]);
    __m128d t2 = _mm_set_pd(a[idx[3]], a[idx[2]]);
    __m128d t3 = _mm_set_pd(a[idx[5]], a[idx[4]]);
    __m128d t4 = _mm_set_pd(a[idx[7]], a[idx[6]]);
    v1 = _mm256_insertf128_pd(v1, t1, 0);
    v1 = _mm256_insertf128_pd(v1, t2, 1);
    v2 = _mm256_insertf128_pd(v2, t3, 0);
    v2 = _mm256_insertf128_pd(v2, t4, 1);
  }

  void scatter(ValueType a[]) const {
    _mm256_storeu_pd(a, v1);
    _mm256_storeu_pd(a+4, v2);
  }

  void aligned_scatter(ValueType a[]) const {
    _mm256_store_pd(a, v1);
    _mm256_store_pd(a+4, v2);
  }

  void zero() {
    v1 = _mm256_setzero_pd();
    v2 = _mm256_setzero_pd();
  }

  void plus_equal(const TinyVec& t) {
    v1 = _mm256_add_pd(v1, t.v1);
    v2 = _mm256_add_pd(v2, t.v2);
  }

  void times_equal(const TinyVec& t) {
    v1 = _mm256_mul_pd(v1, t.v1);
    v2 = _mm256_mul_pd(v2, t.v2);
  }

  // *this = *this + t1 * t2
  void multiply_add(const TinyVec& t1, const TinyVec& t2) {
    __m256d t = _mm256_mul_pd(t1.v1, t2.v1);
    __m256d s = _mm256_mul_pd(t1.v2, t2.v2);
    v1 = _mm256_add_pd(v1, t);
    v2 = _mm256_add_pd(v2, s);
  }

  ValueType sum() {
    __m256d s1 = _mm256_hadd_pd(v1,v1);//v[0]+v[1] v[0]+v[1] v[2]+v[3] v[2]+v[3]
    __m128d s1l = _mm256_extractf128_pd(s1, 0); //v[0]+v[1] v[0]+v[1]
    __m128d s1h = _mm256_extractf128_pd(s1, 1); //v[2]+v[3] v[2]+v[3]
    s1l = _mm_add_pd(s1l,s1h); // v[0]+v[1]+v[2]+v[3] v[0]+v[1]+v[2]+v[3]
    ValueType res1;
    _mm_storel_pd(&res1, s1l);

    __m256d s2 = _mm256_hadd_pd(v2,v2);//v[0]+v[1] v[0]+v[1] v[2]+v[3] v[2]+v[3]
    __m128d s2l = _mm256_extractf128_pd(s2, 0); //v[0]+v[1] v[0]+v[1]
    __m128d s2h = _mm256_extractf128_pd(s2, 1); //v[2]+v[3] v[2]+v[3]
    s2l = _mm_add_pd(s2l,s2h); // v[0]+v[1]+v[2]+v[3] v[0]+v[1]+v[2]+v[3]
    ValueType res2;
    _mm_storel_pd(&res2, s2l);

    return res1 + res2;
  }

private:
  __m256d v1, v2;
};
#endif

#if defined( __MIC__ )
template <>
class TinyVec<double,8,true,false> {
public:

  typedef double ValueType;
  static const int Num = 8;

  TinyVec() {}

  TinyVec(const ValueType a[]) {
    load(a);
  }

  template <typename OrdinalType>
  TinyVec(const ValueType a[], const OrdinalType idx[]) {
    gather(a,idx);
  }

  TinyVec(const ValueType a) {
    load(a);
  }

  TinyVec(const TinyVec& tv) {
    v = tv.v;
  }

  TinyVec& operator=(const TinyVec& tv) {
    v = tv.v;
    return *this;
  }

  void load(const ValueType a[]) {
    v = _mm512_load_pd(a);
  }

  void load(const ValueType a) {
    v = _mm512_set1_pd(a);
  }

  void aligned_load(const ValueType a[]) {
    v = _mm512_load_pd(a);
  }

  template <typename OrdinalType>
  void gather(const ValueType a[], const OrdinalType idx[]) {
    __mmask16 mask = _mm512_int2mask(255);
    __m512i vidx = _mm512_setzero_epi32();
    vidx = _mm512_mask_load_epi32(vidx, mask, idx);
    v = _mm512_i32logather_pd(vidx, a, 8);
  }

  void scatter(ValueType a[]) const {
    _mm512_store_pd(a, v);
  }

  void aligned_scatter(ValueType a[]) const {
    _mm512_store_pd(a, v);
  }

  void zero() {
    v = _mm512_setzero_pd();
  }

  void plus_equal(const TinyVec& t) {
    v = _mm512_add_pd(v, t.v);
  }

  void times_equal(const TinyVec& t) {
    v = _mm512_mul_pd(v, t.v);
  }

  // *this = *this + t1 * t2
  void multiply_add(const TinyVec& t1, const TinyVec& t2) {
    v = _mm512_fmadd_pd(t1.v, t2.v, v);
  }

  ValueType sum() {
    return _mm512_reduce_add_pd(v);
  }

private:
  __m512d v;
};

template <>
class TinyVec<double,8,true,true> {
public:

  typedef double ValueType;
  static const int Num = 8;

  TinyVec(const int sz) {
    mask = _mm512_int2mask((1 << (sz+1))-1);
  }

  TinyVec(const ValueType a[], const int sz) {
    mask = _mm512_int2mask((1 << (sz+1))-1);
    load(a);
  }

  template <typename OrdinalType>
  TinyVec(const ValueType a[], const OrdinalType idx[], const int sz) {
    mask = _mm512_int2mask((1 << (sz+1))-1);
    gather(a,idx);
  }

  TinyVec(const ValueType a, int sz) {
    mask = _mm512_int2mask((1 << (sz+1))-1);
    load(a);
  }

  TinyVec(const TinyVec& tv) {
    mask = tv.mask;
    v = tv.v;
  }

  TinyVec& operator=(const TinyVec& tv) {
    mask = tv.mask;
    v = tv.v;
    return *this;
  }

  void load(const ValueType a[]) {
    v = _mm512_setzero_pd();
    v = _mm512_mask_load_pd(v, mask, a);
  }

  void load(const ValueType a) {
    v = _mm512_set1_pd(a);
  }

  void aligned_load(const ValueType a[]) {
    v = _mm512_setzero_pd();
    v = _mm512_mask_load_pd(v, mask, a);
  }

  template <typename OrdinalType>
  void gather(const ValueType a[], const OrdinalType idx[]) {
    // We're assuming idx is an array of 32-bit integers
    // Load 16 integers into v1idx, then permute the high 256 bits
    // to the low 256 bits (DCBA -> BADC where 128 bit lanes are read right to
    // left).  Then load the vectors into v1 and v2.
    // logather_pd only uses the low 256 bits in the index vector.
    __m512i vidx = _mm512_load_epi32(idx);
    v = _mm512_setzero_pd();
    v = _mm512_mask_i32logather_pd(v, mask, vidx, a, 8);
  }

  void scatter(ValueType a[]) const {
    _mm512_mask_store_pd(a, mask, v);
  }

  void aligned_scatter(ValueType a[]) const {
    _mm512_mask_store_pd(a, mask, v);
  }

  void zero() {
    v = _mm512_setzero_pd();
  }

  void plus_equal(const TinyVec& t) {
    v = _mm512_mask_add_pd(v, mask, v, t.v);
  }

  void times_equal(const TinyVec& t) {
    v = _mm512_mask_mul_pd(v, mask, v, t.v);
  }

  // *this = *this + t1 * t2
  void multiply_add(const TinyVec& t1, const TinyVec& t2) {
    v = _mm512_mask3_fmadd_pd(t1.v, t2.v, v, mask);
  }

  ValueType sum() {
    return _mm512_mask_reduce_add_pd(mask, v);
  }

private:
  __mmask8 mask;
  __m512d v;
};

template <>
class TinyVec<double,16,true,false> {
public:

  typedef double ValueType;
  static const int Num = 16;

  TinyVec() {}

  TinyVec(const ValueType a[]) {
    load(a);
  }

  template <typename OrdinalType>
  TinyVec(const ValueType a[], const OrdinalType idx[]) {
    gather(a,idx);
  }

  TinyVec(const ValueType a) {
    load(a);
  }

  TinyVec(const TinyVec& tv) {
    v1 = tv.v1; v2 = tv.v2;
  }

  TinyVec& operator=(const TinyVec& tv) {
    v1 = tv.v1; v2 = tv.v2;
    return *this;
  }

  void load(const ValueType a[]) {
    v1 = _mm512_load_pd(a);
    v2 = _mm512_load_pd(a+8);
  }

  void load(const ValueType a) {
    v1 = _mm512_set1_pd(a);
    v2 = _mm512_set1_pd(a);
  }

  void aligned_load(const ValueType a[]) {
    v1 = _mm512_load_pd(a);
    v2 = _mm512_load_pd(a+8);
  }

  template <typename OrdinalType>
  void gather(const ValueType a[], const OrdinalType idx[]) {
    // We're assuming idx is an array of 32-bit integers
    // Load 16 integers into v1idx, then permute the high 256 bits
    // to the low 256 bits (DCBA -> BADC where 128 bit lanes are read right to
    // left).  Then load the vectors into v1 and v2.
    // logather_pd only uses the low 256 bits in the index vector.
    __m512i v1idx = _mm512_load_epi32(idx);
    __m512i v2idx = _mm512_permute4f128_epi32(v1idx, _MM_PERM_BADC);
    v1 = _mm512_i32logather_pd(v1idx, a, 8);
    v2 = _mm512_i32logather_pd(v2idx, a, 8);
  }

  void scatter(ValueType a[]) const {
    _mm512_store_pd(a, v1);
    _mm512_store_pd(a+8, v2);
  }

  void aligned_scatter(ValueType a[]) const {
    _mm512_store_pd(a, v1);
    _mm512_store_pd(a+8, v2);
  }

  void zero() {
    v1 = _mm512_setzero_pd();
    v2 = _mm512_setzero_pd();
  }

  void plus_equal(const TinyVec& t) {
    v1 = _mm512_add_pd(v1, t.v1);
    v2 = _mm512_add_pd(v2, t.v2);
  }

  void times_equal(const TinyVec& t) {
    v1 = _mm512_mul_pd(v1, t.v1);
    v2 = _mm512_mul_pd(v2, t.v2);
  }

  // *this = *this + t1 * t2
  void multiply_add(const TinyVec& t1, const TinyVec& t2) {
    v1 = _mm512_fmadd_pd(t1.v1, t2.v1, v1);
    v2 = _mm512_fmadd_pd(t1.v2, t2.v2, v2);
  }

  ValueType sum() {
    return _mm512_reduce_add_pd(v1) + _mm512_reduce_add_pd(v2);
  }

private:
  __m512d v1, v2;
};

template <>
class TinyVec<double,16,true,true> {
public:

  typedef double ValueType;
  static const int Num = 16;

  TinyVec(const int sz) {
    mask = _mm512_int2mask((1 << (sz-7))-1);
  }

  TinyVec(const ValueType a[], int sz) {
    mask = _mm512_int2mask((1 << (sz-7))-1);
    load(a);
  }

  template <typename OrdinalType>
  TinyVec(const ValueType a[], const OrdinalType idx[], int sz) {
    mask = _mm512_int2mask((1 << (sz-7))-1);
    gather(a,idx);
  }

  TinyVec(const ValueType a, int sz) {
    mask = _mm512_int2mask((1 << (sz-7))-1);
    load(a);
  }

  TinyVec(const TinyVec& tv) {
    mask = tv.mask;
    v1 = tv.v1; v2 = tv.v2;
  }

  TinyVec& operator=(const TinyVec& tv) {
    mask = tv.mask;
    v1 = tv.v1; v2 = tv.v2;
    return *this;
  }

  void load(const ValueType a[]) {
    v1 = _mm512_load_pd(a);
    v2 = _mm512_setzero_pd();
    v2 = _mm512_mask_load_pd(v2, mask, a+8);
  }

  void load(const ValueType a) {
    v1 = _mm512_set1_pd(a);
    v2 = _mm512_set1_pd(a);
  }

  void aligned_load(const ValueType a[]) {
    v1 = _mm512_load_pd(a);
    v2 = _mm512_setzero_pd();
    v2 = _mm512_mask_load_pd(v2, mask, a+8);
  }

  template <typename OrdinalType>
  void gather(const ValueType a[], const OrdinalType idx[]) {
    // We're assuming idx is an array of 32-bit integers
    // Load 16 integers into v1idx, then permute the high 256 bits
    // to the low 256 bits (DCBA -> BADC where 128 bit lanes are read right to
    // left).  Then load the vectors into v1 and v2.
    // logather_pd only uses the low 256 bits in the index vector.
    // Note:  permute4f128 overwrites its argument, so we need to load v1 first
    __m512i v1idx = _mm512_load_epi32(idx);
    v1 = _mm512_i32logather_pd(v1idx, a, 8);

    v1idx = _mm512_permute4f128_epi32(v1idx, _MM_PERM_BADC);
    v2 = _mm512_setzero_pd();
    v2 = _mm512_mask_i32logather_pd(v2, mask, v1idx, a, 8);
  }

  void scatter(ValueType a[]) const {
    _mm512_store_pd(a, v1);
    _mm512_mask_store_pd(a+8, mask, v2);
  }

  void aligned_scatter(ValueType a[]) const {
    _mm512_store_pd(a, v1);
    _mm512_mask_store_pd(a+8, mask, v2);
  }

  void zero() {
    v1 = _mm512_setzero_pd();
    v2 = _mm512_setzero_pd();
  }

  void plus_equal(const TinyVec& t) {
    v1 = _mm512_add_pd(v1, t.v1);
    v2 = _mm512_mask_add_pd(v2, mask, v2, t.v2);
  }

  void times_equal(const TinyVec& t) {
    v1 = _mm512_mul_pd(v1, t.v1);
    v2 = _mm512_mask_mul_pd(v2, mask, v2, t.v2);
  }

  // *this = *this + t1 * t2
  void multiply_add(const TinyVec& t1, const TinyVec& t2) {
    v1 = _mm512_fmadd_pd(t1.v1, t2.v1, v1);
    v2 = _mm512_mask3_fmadd_pd(t1.v2, t2.v2, v2, mask);
  }

  ValueType sum() {
    return _mm512_reduce_add_pd(v1) + _mm512_mask_reduce_add_pd(mask, v2);
  }

private:
  __mmask8 mask;
  __m512d v1, v2;
};
#endif

#endif // #if defined(HAVE_STOKHOS_INTRINSICS) && !defined( __CUDACC__ )

} // namespace Stokhos

#endif /* #ifndef STOKHOS_TINY_VEC_HPP */

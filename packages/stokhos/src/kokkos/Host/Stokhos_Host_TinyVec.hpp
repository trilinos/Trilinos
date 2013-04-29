// @HEADER
// ***********************************************************************
//
//                     Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_HOST_TINY_VEC_HPP
#define STOKHOS_HOST_TINY_VEC_HPP

extern "C" {
#include <immintrin.h>
}

namespace Stokhos {

template <typename ValueType, int N, bool UseIntrinsics >
class TinyVec {
public:

  static const int Num = N;

  TinyVec() {}

  template <typename ValueArrayType>
  TinyVec(const ValueArrayType a) {
    load(a);
  }

  template <typename ValueArrayType, typename OrdinalArrayType>
  TinyVec(const ValueArrayType a, const OrdinalArrayType idx) {
    gather(a,idx);
  }

  template <typename ValueArrayType, typename OrdinalArrayType,
            typename OrdinalType>
  TinyVec(const ValueArrayType a, const OrdinalArrayType idx,
          const OrdinalType stride) {
    gather(a,idx,stride);
  }

  TinyVec(const TinyVec& tv) {
    for (int i=0; i<Num; ++i)
      v[i] = tv.v[i];
  }

  TinyVec& operator=(const TinyVec& tv) {
    for (int i=0; i<Num; ++i)
      v[i] = tv.v[i];
    return *this;
  }

  template <typename ValueArrayType>
  void load(const ValueArrayType a) {
    for (int i=0; i<Num; ++i)
      v[i] = a[i];
  }

  void load(const ValueType a) {
    for (int i=0; i<Num; ++i)
      v[i] = a;
  }

  template <typename ValueArrayType>
  void aligned_load(const ValueArrayType a) {
    for (int i=0; i<Num; ++i)
      v[i] = a[i];
  }

  template <typename ValueArrayType, typename OrdinalArrayType>
  void gather(const ValueArrayType a, const OrdinalArrayType idx) {
    for (int i=0; i<Num; ++i)
      v[i] = a[idx[i]];
  }

  template <typename ValueArrayType, typename OrdinalArrayType,
            typename OrdinalType>
  void gather(const ValueArrayType a, const OrdinalArrayType idx,
              const OrdinalType stride) {
    for (int i=0; i<Num; ++i)
      v[i] = a[idx[i*stride]];
  }

  template <typename ValueArrayType>
  void scatter(ValueArrayType a) const {
    for (int i=0; i<Num; ++i)
      a[i] = v[i];
  }

  template <typename ValueArrayType>
  void aligned_scatter(ValueArrayType a) const {
    for (int i=0; i<Num; ++i)
      a[i] = v[i];
  }

  void zero() {
    for (int i=0; i<Num; ++i)
      v[i] = ValueType(0.0);
  }

  void plus_equal(const TinyVec& t) {
    for (int i=0; i<Num; ++i)
      v[i] += t.v[i];
  }

  void times_equal(const TinyVec& t) {
    for (int i=0; i<Num; ++i)
      v[i] *= t.v[i];
  }

  ValueType sum() const {
    ValueType s(0.0);
    for (int i=0; i<Num; ++i)
      s += v[i];
    return s;
  }

private:
   ValueType v[Num];
};

#ifdef __SSE2__
template <>
class TinyVec<double,2,true> {
public:

  typedef double ValueType;
  static const int Num = 2;

  TinyVec() {}

  template <typename ValueArrayType>
  TinyVec(const ValueArrayType a) {
    load(a);
  }

  template <typename ValueArrayType, typename OrdinalArrayType>
  TinyVec(const ValueArrayType a, const OrdinalArrayType idx) {
    gather(a,idx);
  }

  template <typename ValueArrayType, typename OrdinalArrayType,
            typename OrdinalType>
  TinyVec(const ValueArrayType a, const OrdinalArrayType idx,
          const OrdinalType stride) {
    gather(a,idx,stride);
  }

  TinyVec(const TinyVec& tv) {
    v = tv.v;
  }

  TinyVec& operator=(const TinyVec& tv) {
    v = tv.v;
    return *this;
  }

  template <typename ValueArrayType>
  void load(const ValueArrayType a) {
    v = _mm_set_pd(a[1], a[0]);
  }

  void load(const ValueType a) {
    v = _mm_set1_pd(a);
  }

  template <typename ValueArrayType>
  void aligned_load(const ValueArrayType a) {
    v = _mm_load_pd(a);
  }

  template <typename ValueArrayType, typename OrdinalArrayType>
  void gather(const ValueArrayType a, const OrdinalArrayType idx) {
    v = _mm_set_pd(a[idx[1]], a[idx[0]]);
  }

  template <typename ValueArrayType, typename OrdinalArrayType,
            typename OrdinalType>
  void gather(const ValueArrayType a, const OrdinalArrayType idx,
              const OrdinalType stride) {
    v = _mm_set_pd(a[idx[stride]], a[idx[0]]);
  }

  template <typename ValueArrayType>
  void scatter(ValueArrayType a) const {
    _mm_storel_pd(&a[0], v);
    _mm_storeh_pd(&a[1], v);
  }

  template <typename ValueArrayType>
  void aligned_scatter(ValueArrayType a) const {
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

  ValueType sum() const {
    ValueType a[Num];
    aligned_scatter(a);
    return a[0]+a[1];
  }

private:
  __m128d v;
};
#endif

//#ifdef __AVX__
template <>
class TinyVec<double,4,true> {
public:

  typedef double ValueType;
  static const int Num = 4;

  TinyVec() {}

  template <typename ValueArrayType>
  TinyVec(const ValueArrayType a) {
    load(a);
  }

  template <typename ValueArrayType, typename OrdinalArrayType>
  TinyVec(const ValueArrayType a, const OrdinalArrayType idx) {
    gather(a,idx);
  }

  template <typename ValueArrayType, typename OrdinalArrayType,
            typename OrdinalType>
  TinyVec(const ValueArrayType a, const OrdinalArrayType idx,
          const OrdinalType stride) {
    gather(a,idx,stride);
  }

  TinyVec(const TinyVec& tv) {
    v = tv.v;
  }

  TinyVec& operator=(const TinyVec& tv) {
    v = tv.v;
    return *this;
  }

  template <typename ValueArrayType>
  void load(const ValueArrayType a) {
    v = _mm256_loadu_pd(a);
  }

  void load(const ValueType a) {
    v = _mm256_set1_pd(a);
  }

  template <typename ValueArrayType>
  void aligned_load(const ValueArrayType a) {
    v = _mm256_load_pd(a);
  }

  template <typename ValueArrayType, typename OrdinalArrayType>
  void gather(const ValueArrayType a, const OrdinalArrayType idx) {
    __m128d v1 = _mm_set_pd(a[idx[1]], a[idx[0]]);
    __m128d v2 = _mm_set_pd(a[idx[3]], a[idx[2]]);
    v = _mm256_insertf128_pd(v, v1, 0);
    v = _mm256_insertf128_pd(v, v2, 1);
  }

  template <typename ValueArrayType, typename OrdinalArrayType,
            typename OrdinalType>
  void gather(const ValueArrayType a, const OrdinalArrayType idx,
              const OrdinalType stride) {
    __m128d v1 = _mm_set_pd(a[idx[stride]], a[idx[0]]);
    __m128d v2 = _mm_set_pd(a[idx[3*stride]], a[idx[2*stride]]);
    v = _mm256_insertf128_pd(v, v1, 0);
    v = _mm256_insertf128_pd(v, v2, 1);
  }

  template <typename ValueArrayType>
  void scatter(ValueArrayType a) {
    _mm256_storeu_pd(a, v);
  }

  template <typename ValueArrayType>
  void aligned_scatter(ValueArrayType a) {
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
//#endif

} // namespace Stokhos

#endif /* #ifndef STOKHOS_HOST_TINY_VEC_HPP */

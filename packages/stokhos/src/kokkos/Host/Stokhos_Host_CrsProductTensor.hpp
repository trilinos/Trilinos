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

#ifndef STOKHOS_HOST_CRSPRODUCTTENSOR_HPP
#define STOKHOS_HOST_CRSPRODUCTTENSOR_HPP

#include "KokkosArray_Host.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_CrsProductTensor.hpp"

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

  template <typename ValueArrayType>
  void load(const ValueArrayType a) {
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

  template <typename ValueArrayType>
  void load(const ValueArrayType a) {
    //v = _mm_load_pd(a);
    v = _mm_set_pd(a[0], a[1]);
  }

  template <typename ValueArrayType, typename OrdinalArrayType>
  void gather(const ValueArrayType a, const OrdinalArrayType idx) {
    v = _mm_set_pd(a[idx[0]], a[idx[1]]);
  }

  template <typename ValueArrayType, typename OrdinalArrayType,
	    typename OrdinalType>
  void gather(const ValueArrayType a, const OrdinalArrayType idx,
	      const OrdinalType stride) {
    v = _mm_set_pd(a[idx[0]], a[idx[stride]]);
  }

  template <typename ValueArrayType>
  void scatter(ValueArrayType a) const {
    //_mm_store_pd(a, v);
    _mm_storel_pd(&a[0], v);
    _mm_storeh_pd(&a[1], v);
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
    scatter(a);
    return a[0]+a[1];
  }

private:
  __m128d v;
};
#endif

#ifdef __AVX__
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

  template <typename ValueArrayType>
  void load(const ValueArrayType a) {
    v = _mm256_load_pd(a);
  }

  template <typename ValueArrayType, typename OrdinalArrayType>
  void gather(const ValueArrayType a, const OrdinalArrayType idx) {
    __m128d v1 = _mm_set_pd(a[idx[0]], a[idx[1]]);
    __m128d v2 = _mm_set_pd(a[idx[2]], a[idx[3]]);
    v = _mm256_insertf128_pd(v, v1, 0);
    v = _mm256_insertf128_pd(v, v2, 1);
  }

  template <typename ValueArrayType, typename OrdinalArrayType,
	    typename OrdinalType>
  void gather(const ValueArrayType a, const OrdinalArrayType idx,
	      const OrdinalType stride) {
    __m128d v1 = _mm_set_pd(a[idx[0]], a[idx[stride]]);
    __m128d v2 = _mm_set_pd(a[idx[2*stride]], a[idx[3*stride]]);
    v = _mm256_insertf128_pd(v, v1, 0);
    v = _mm256_insertf128_pd(v, v2, 1);
  }

   template <typename ValueArrayType>
  void scatter(ValueArrayType a) {
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
    ValueType a[Num];
    scatter(a);
    return a[0]+a[1]+a[2]+a[3];
  }

private:
  __m256d v;
};
#endif

template< typename ValueType >
class Multiply< CrsProductTensor< ValueType , KokkosArray::Host > , void , void , DefaultSparseMatOps >
{
public:
  
  typedef KokkosArray::Host::size_type size_type ;
  typedef CrsProductTensor< ValueType , KokkosArray::Host > tensor_type ;

  template< typename MatrixValue , typename VectorValue >
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {
    const size_type block_size = 2;
    typedef TinyVec<ValueType,block_size,false> TV;

    const size_type nDim = tensor.dimension();    

    for ( size_type iy = 0 ; iy < nDim ; ++iy ) {

      const size_type nEntry = tensor.num_entry(iy);
      const size_type iEntryBeg = tensor.entry_begin(iy);
      const size_type iEntryEnd = iEntryBeg + nEntry;
            size_type iEntry    = iEntryBeg;

      VectorValue ytmp = 0 ;

      // Do entries with a blocked loop of size block_size
      if (block_size > 1) {
	const size_type nBlock = nEntry / block_size;
	const size_type nEntryB = nBlock * block_size;
	const size_type iEnd = iEntryBeg + nEntryB;
	
	TV vy;
	vy.zero();
	int j[block_size], k[block_size];
	
	for ( ; iEntry < iEnd ; iEntry += block_size ) {
	  
	  for (size_type ii=0; ii<block_size; ++ii) {
	    j[ii] = tensor.coord(iEntry+ii,0);
	    k[ii] = tensor.coord(iEntry+ii,1);
	  }
	  TV aj(a, j), ak(a, k), xj(x, j), xk(x, k), 
	    c(&(tensor.value(iEntry)));
	  
	  // const size_type *j = &(tensor.coord(iEntry,0));
	  // const size_type *k = &(tensor.coord(iEntry,1));
	  // TV aj(a, j, 2), ak(a, k, 2), xj(x, j, 2), xk(x, k, 2), 
	  //    c(&(tensor.value(iEntry)));
	  
	  // vy += c * ( aj * xk + ak * xj)
	  aj.times_equal(xk);
	  ak.times_equal(xj);
	  aj.plus_equal(ak);
	  c.times_equal(aj);
	  vy.plus_equal(c);
	  
	}
	
	ytmp += vy.sum();
      }
      
      // Do remaining entries with a scalar loop
      for ( ; iEntry<iEntryEnd; ++iEntry) {
	const size_type j = tensor.coord(iEntry,0);
	const size_type k = tensor.coord(iEntry,1);
        
	ytmp += tensor.value(iEntry) * ( a[j] * x[k] + a[k] * x[j] );
      }

      y[iy] += ytmp ;
    }
  }

  static size_type matrix_size( const tensor_type & tensor )
  { return tensor.dimension(); }

  static size_type vector_size( const tensor_type & tensor )
  { return tensor.dimension(); }
};

//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_HOST_SPARSEPRODUCTTENSOR_HPP */


/*
//@HEADER
// ************************************************************************
// 
//    KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_HOST_SPARSEPRODUCTTENSOR_HPP
#define KOKKOSARRAY_HOST_SPARSEPRODUCTTENSOR_HPP

extern "C" {
#include <immintrin.h>
}

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< typename ValueType >
class Multiply< SparseProductTensor< 3 , ValueType , Host > , void , void >
{
public:

  typedef Host::size_type size_type ;
  typedef SparseProductTensor< 3 , ValueType , Host > tensor_type ;

  template< typename MatrixValue , typename VectorValue >
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {
    const size_type nEntry = tensor.entry_count();

    for ( size_type iEntry = 0 ; iEntry < nEntry ; ++iEntry ) {
      const size_type i = tensor.coord(iEntry,0);
      const size_type j = tensor.coord(iEntry,1);
      const size_type k = tensor.coord(iEntry,2);
      const ValueType v = tensor.value(iEntry);

      const bool neq_ij = i != j ;
      const bool neq_jk = j != k ;
      const bool neq_ki = k != i ;

      y[k] += neq_ij ? v * ( a[i] * x[j] + x[i] * a[j] )
                     : v * ( a[i] * x[i] );

      if ( neq_jk ) {
        y[j] += neq_ki ? v * ( a[i] * x[k] + x[i] * a[k] )
                       : v * ( a[i] * x[i] );
      }

      if ( neq_ki && neq_ij ) {
        y[i] += neq_jk ? v * ( a[k] * x[j] + x[k] * a[j] )
                       : v * ( a[j] * x[j] );
      }
    }
  }

  static size_type matrix_size( const tensor_type & tensor )
  { return tensor.dimension(); }

  static size_type vector_size( const tensor_type & tensor )
  { return tensor.dimension(); }
};

//----------------------------------------------------------------------------

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
class Multiply< CrsProductTensor< 3 , ValueType , Host > , void , void >
{
public:

  typedef Host::size_type size_type ;
  typedef CrsProductTensor< 3 , ValueType , Host > tensor_type ;

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

template< typename ValueType >
class Multiply< SparseProductTensor< 3 , ValueType , Host > ,
                SymmetricDiagonalSpec< KokkosArray::Host > , void >
{
public:
  typedef Host::size_type size_type ;
  typedef SparseProductTensor< 3 , ValueType , Host > tensor_type ;

  template< typename MatrixValue >
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                           MatrixValue * const M )
  {
    const SymmetricDiagonalSpec< KokkosArray::Host >
      spec( tensor.dimension() );

    const size_type nEntry  = tensor.entry_count();
    const size_type nMatrix = spec.matrix_size();

    for ( size_type iMatrix = 0 ; iMatrix < nMatrix ; ++iMatrix ) {
      M[ iMatrix ] = 0 ;
    }

    for ( size_type iEntry = 0 ; iEntry < nEntry ; ++iEntry ) {
      const size_type i = tensor.coord(iEntry,0);
      const size_type j = tensor.coord(iEntry,1);
      const size_type k = tensor.coord(iEntry,2);
      const ValueType v = tensor.value(iEntry);

      M[ spec.matrix_offset(k,j) ] += v * a[i] ;

      if ( i != j ) {
        M[ spec.matrix_offset(k,i) ] += v * a[j] ;
      }

      if ( i != k && j != k ) {
        M[ spec.matrix_offset(j,i) ] += v * a[k] ;
      }
    }
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_HOST_SPARSEPRODUCTTENSOR_HPP */


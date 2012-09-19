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

#ifndef KOKKOSARRAY_PRODUCTTENSOR_HPP
#define KOKKOSARRAY_PRODUCTTENSOR_HPP

#include <map>
#include <KokkosArray_CrsArray.hpp>
#include <impl/KokkosArray_Multiply.hpp>

namespace KokkosArray {

//----------------------------------------------------------------------------
/** \brief  Use symmetry to compress the product tensor index space.
 *
 *  By symmetry the coordinates of a particular index can be
 *  arbitrarily reordered.  Compression due to symmetry results
 *  in { D! / ( (D-Rank)! * Rank! ) } unique entries in the index space.
 *
 *  Indices are lexigraphically ordered.
 *  coord(0) >= coord(1) >= ... >= coord(Rank-1)
 */
template< unsigned Rank , class Device >
class ProductTensorIndex {
public:
  typedef Device                           device_type ;
  typedef typename device_type::size_type  size_type ;

  ProductTensorIndex();
  ProductTensorIndex( const ProductTensorIndex & );
  ProductTensorIndex & operator = ( const ProductTensorIndex & );

  explicit ProductTensorIndex( size_type offset );

  // ProductTensorIndex( size_type coord0 , size_type coord1 , ... );

  ProductTensorIndex & operator ++ ();

  /** \brief  Coordinate 'c' where 0 <= c < Rank */
  size_type coord( size_type c ) const ;

  /** \brief  Offset of this entry in the index space */
  size_type offset() const ;

  bool operator == ( const ProductTensorIndex & ) const ;
  bool operator != ( const ProductTensorIndex & ) const ;
  bool operator <  ( const ProductTensorIndex & ) const ;
  bool operator <= ( const ProductTensorIndex & ) const ;
  bool operator >  ( const ProductTensorIndex & ) const ;
  bool operator >= ( const ProductTensorIndex & ) const ;
};

//----------------------------------------------------------------------------

/** \brief  Sparse product tensor using coordinate storage.
 */
template< unsigned Rank , typename ValueType , class Device >
class SparseProductTensor {
public:
  typedef Device                           device_type ;
  typedef typename device_type::size_type  size_type ;
  typedef ValueType                        value_type ;

  ~SparseProductTensor();
  SparseProductTensor();
  SparseProductTensor( const SparseProductTensor & );
  SparseProductTensor & operator = ( const SparseProductTensor & );

  /** \brief  Dimension of the tensor. */
  size_type dimension() const ;

  /** \brief  Number of sparse entries. */
  size_type entry_count() const ;

  /** \brief  Coordinates of an entry */
  size_type coord( size_type offset , size_type c ) const ;

  /** \brief  Value of an entry */
  const value_type & value( size_type offset ) const ;

private:

  template< class Tensor , class Input >
  friend class CreateSparseProductTensor ;
};

// Specialization for Multiply action of a tensor:
//
// Multiply< SparseProductTensor<R,V,D> >
//   ::apply( const SparseProductTensor<R,V,D> & block ,
//            const MatrixValue       * A ,
//            const InputVectorValue  * x ,
//                  OutputVectorValue * y );
// Where:
//
//  A[ Multiply< SparseProductTensor<R,V,D> >::matrix_size( block ) ]
//  x[ Multiply< SparseProductTensor<R,V,> >::vector_size( block ) ]
//  y[ Multiply< SparseProductTensor<R,V,> >::vector_size( block ) ]
//----------------------------------------------------------------------------

/** \brief  Sparse product tensor with replicated entries
 *          to provide subsets with a given coordinate.
 *
 *  This allows product tensor multiplication to be partitioned
 *  on a given coordinate values.
 *
 *  for ( size_type i = 0 ; i < p.dimension() ; ++i ) {
 *    y[i] = 0 ;
 *    for ( size_type e = p.entry_begin(i) ; 
 *                    e < p.entry_end(i) ; ++e ) {
 *      const size_type j = p.coord(e,0);
 *      const size_type k = p.coord(e,1);
 *      Scalar tmp = a[j] * x[k] ; if ( j != k ) tmp += a[k] * x[j] ;
 *      y[i] += p.value(e) * tmp ;
 *    }
 *  }
 */
template< unsigned Rank , typename ValueType , class Device >
class CrsProductTensor {
public:
  typedef Device                           device_type ;
  typedef typename device_type::size_type  size_type ;
  typedef ValueType                        value_type ;

  ~CrsProductTensor();
  CrsProductTensor();
  CrsProductTensor( const CrsProductTensor & );
  CrsProductTensor & operator = ( const CrsProductTensor & );

  /** \brief  Dimension of the tensor. */
  size_type dimension() const ;

  /** \brief  Number of sparse entries. */
  size_type entry_count() const ;

  /** \brief  Maximum sparse entries for any coordinate */
  size_type entry_maximum() const ;

  /** \brief  Begin entries with a coordinate 'i' */
  size_type entry_begin( size_type i ) const ;

  /** \brief  End entries with a coordinate 'i' */
  size_type entry_end( size_type i ) const ;

  /** \brief  Coordinates of an entry */
  size_type coord( size_type entry , size_type c ) const ;

  /** \brief  Value of an entry */
  const value_type & value( size_type entry ) const ;
};

//----------------------------------------------------------------------------

namespace Impl {

template< class Tensor , class Input >
class CreateProductTensor ;

template< class Tensor , class Input >
class CreateSparseProductTensor ;

}

template< class Tensor , class Input >
typename Impl::CreateSparseProductTensor<Tensor,Input>::type
create_product_tensor( const Input & input )
{
  typedef typename Impl::CreateSparseProductTensor<Tensor,Input> factory ;
  return factory::create( input );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace KokkosArray

#include <impl/KokkosArray_ProductTensor_create.hpp>

#endif /* #ifndef KOKKOSARRAY_PRODUCTTENSOR_HPP */



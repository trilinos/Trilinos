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

#include <KokkosArray_ProductTensorIndex.hpp>
#include <KokkosArray_CrsArray.hpp>
#include <impl/KokkosArray_Multiply.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class Tensor , class Input >
class CreateProductTensor ;

template< class Tensor , class Input >
class CreateSparseProductTensor ;

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

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

template< typename ValueType , class DeviceType >
class SparseProductTensor< 3 , ValueType , DeviceType > {
public:

  typedef DeviceType                       device_type ;
  typedef typename device_type::size_type  size_type ;
  typedef ValueType                        value_type ;

private:

  typedef KokkosArray::View< size_type[][3] , device_type >  map_type ;
  typedef KokkosArray::View< value_type[] ,   device_type >  vec_type ;

public:

  inline
  ~SparseProductTensor() {}

  inline
  SparseProductTensor() : m_coord() , m_value() , m_dimen(0) {}

  inline
  SparseProductTensor( const SparseProductTensor & rhs )
  : m_coord( rhs.m_coord ) , m_value( rhs.m_value ) , m_dimen( rhs.m_dimen ) {}

  inline
  SparseProductTensor & operator = ( const SparseProductTensor & rhs )
  {
    m_coord = rhs.m_coord ;
    m_value = rhs.m_value ;
    m_dimen = rhs.m_dimen ;
    return *this ;
  }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension() const { return m_dimen ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type entry_count() const
  { return m_value.dimension_0(); }

  KOKKOSARRAY_INLINE_FUNCTION
  const value_type & value( const size_type entry ) const
  { return m_value( entry ); }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type coord( const size_type entry , const size_type c ) const
  { return m_coord( entry , c ); }

private:

  map_type   m_coord ;
  vec_type   m_value ;
  size_type  m_dimen ;

  template< class T , class I >
  friend class Impl::CreateSparseProductTensor ;
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

} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

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

template< typename ValueType , class DeviceType >
class CrsProductTensor< 3 , ValueType , DeviceType > {
public:

  typedef DeviceType                       device_type ;
  typedef typename device_type::size_type  size_type ;
  typedef ValueType                        value_type ;

private:

  typedef KokkosArray::View< value_type[] , device_type >  vec_type ;

  KokkosArray::CrsArray< size_type[2] , device_type >  m_coord ;
  KokkosArray::View< value_type[] , device_type >      m_value ;
  size_type                                            m_entry_max ;

  template< class T , class I >
  friend class Impl::CreateSparseProductTensor ;

public:

  inline
  ~CrsProductTensor() {}

  inline
  CrsProductTensor() : m_coord() , m_value() , m_entry_max(0) {}

  inline
  CrsProductTensor( const CrsProductTensor & rhs )
  : m_coord( rhs.m_coord ) , m_value( rhs.m_value ) , m_entry_max( rhs.m_entry_max ) {}

  inline
  CrsProductTensor & operator = ( const CrsProductTensor & rhs )
  {
    m_coord = rhs.m_coord ;
    m_value = rhs.m_value ;
    m_entry_max = rhs.m_entry_max ;
    return *this ;
  }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension() const { return m_coord.row_map.dimension(0) - 1 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type entry_count() const
  { return m_coord.entries.dimension(0); }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type entry_maximum() const
  { return m_entry_max ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type entry_begin( size_type i ) const
  { return m_coord.row_map[i]; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type entry_end( size_type i ) const
  { return m_coord.row_map[i+1]; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type coord( const size_type entry , const size_type c ) const
  { return m_coord.entries( entry , c ); }

  KOKKOSARRAY_INLINE_FUNCTION
  const value_type & value( const size_type entry ) const
  { return m_value( entry ); }

};

} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

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



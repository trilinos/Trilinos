// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_FLAT_SPARSE_3_TENSOR_HPP
#define STOKHOS_FLAT_SPARSE_3_TENSOR_HPP

#include "KokkosArray_View.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Stokhos {

/** \brief  Sparse product tensor with replicated entries
 *          to provide subsets with a given coordinate.
 */
template< typename ValueType , class DeviceType >
class FlatSparse3Tensor {
public:

  typedef DeviceType                       device_type ;
  typedef typename device_type::size_type  size_type ;
  typedef ValueType                        value_type ;

private:

  typedef KokkosArray::View< size_type[] , device_type > coord_array_type ;
  typedef KokkosArray::View< value_type[], device_type > value_array_type ;
  typedef KokkosArray::View< size_type[], device_type > entry_array_type ;
  typedef KokkosArray::View< size_type[], device_type > row_map_array_type ;

  coord_array_type   m_k_coord ;
  coord_array_type   m_j_coord ;
  value_array_type   m_value ;
  entry_array_type   m_num_k ;
  entry_array_type   m_num_j ;
  row_map_array_type m_k_row_map ;
  row_map_array_type m_j_row_map ;
  size_type          m_nnz ;
  size_type          m_flops ;

public:

  inline
  ~FlatSparse3Tensor() {}

  inline
  FlatSparse3Tensor() : 
    m_k_coord() , 
    m_j_coord() , 
    m_value() , 
    m_num_k() , 
    m_num_j() , 
    m_k_row_map() , 
    m_j_row_map() , 
    m_nnz(0) ,
    m_flops(0) {}

  inline
  FlatSparse3Tensor( const FlatSparse3Tensor & rhs ) : 
    m_k_coord( rhs.m_k_coord ) , 
    m_j_coord( rhs.m_j_coord ) ,
    m_value( rhs.m_value ) , 
    m_num_k( rhs.m_num_k ) ,
    m_num_j( rhs.m_num_j ) ,
    m_k_row_map( rhs.m_k_row_map ) , 
    m_j_row_map( rhs.m_j_row_map ) , 
    m_nnz( rhs.m_nnz ) ,
    m_flops( rhs.m_flops ) {}

  inline
  FlatSparse3Tensor & operator = ( const FlatSparse3Tensor & rhs )
  {
    m_k_coord = rhs.m_k_coord ;
    m_j_coord = rhs.m_j_coord ;
    m_value = rhs.m_value ;
    m_num_k = rhs.m_num_k ;
    m_num_j = rhs.m_num_j ;
    m_k_row_map = rhs.m_k_row_map ;
    m_j_row_map = rhs.m_j_row_map ;
    m_nnz = rhs.m_nnz;
    m_flops = rhs.m_flops;
    return *this ;
  }

  /** \brief  Dimension of the tensor. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension() const { return m_k_row_map.dimension(0) - 1 ; }

  /** \brief  Number of sparse entries. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type entry_count() const
  { return m_j_coord.dimension(0); }

  /** \brief  Begin k entries with a coordinate 'i' */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type k_begin( size_type i ) const
  { return m_k_row_map[i]; }

  /** \brief  End k entries with a coordinate 'i' */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type k_end( size_type i ) const
  { return m_k_row_map[i] + m_num_k(i); }

  /** \brief  Number of k entries with a coordinate 'i' */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_k( size_type i ) const
  { return m_num_k(i); }

  /** \brief  k coordinate for k entry 'kEntry' */
  KOKKOSARRAY_INLINE_FUNCTION
  const size_type& k_coord( const size_type kEntry ) const
  { return m_k_coord( kEntry ); }

  /** \brief  Begin j entries with a k entry 'kEntry' */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type j_begin( size_type kEntry ) const
  { return m_j_row_map[kEntry]; }

  /** \brief  End j entries with a k entry 'kEntry' */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type j_end( size_type kEntry ) const
  { return m_j_row_map[kEntry] + m_num_j(kEntry); }

  /** \brief  Number of j entries with a k entry 'kEntry' */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_j( size_type kEntry ) const
  { return m_num_j(kEntry); }

  /** \brief  j coordinate for j entry 'jEntry' */
  KOKKOSARRAY_INLINE_FUNCTION
  const size_type& j_coord( const size_type jEntry ) const
  { return m_j_coord( jEntry ); }

  /** \brief  Value for j entry 'jEntry' */
  KOKKOSARRAY_INLINE_FUNCTION
  const value_type & value( const size_type jEntry ) const
  { return m_value( jEntry ); }

  /** \brief Number of non-zero's */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_non_zeros() const 
  { return m_nnz; }

  /** \brief Number flop's per multiply-add */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_flops() const 
  { return m_flops; }

  template <typename OrdinalType>
  static FlatSparse3Tensor
  create( const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
	  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk )
  {
    typedef Stokhos::Sparse3Tensor<OrdinalType,ValueType> Cijk_type;
    
    // Compute number of k's for each i
    const size_type dimension = basis.size();
    std::vector< size_t > k_coord_work( dimension , (size_t) 0 );
    size_type k_entry_count = 0 ;
    for (typename Cijk_type::i_iterator i_it=Cijk.i_begin(); 
	 i_it!=Cijk.i_end(); ++i_it) {
      OrdinalType i = index(i_it);
      k_coord_work[i] = Cijk.num_k(i_it);
      k_entry_count += Cijk.num_k(i_it);
    }

    // Compute number of j's for each i and k
    std::vector< size_t > j_coord_work( k_entry_count , (size_t) 0 );
    size_type j_entry_count = 0 ;
    size_type k_entry = 0 ;
    for (typename Cijk_type::i_iterator i_it=Cijk.i_begin(); 
	 i_it!=Cijk.i_end(); ++i_it) {
      for (typename Cijk_type::ik_iterator k_it = Cijk.k_begin(i_it); 
	   k_it != Cijk.k_end(i_it); ++k_it, ++k_entry) {
	OrdinalType k = index(k_it);
	for (typename Cijk_type::ikj_iterator j_it = Cijk.j_begin(k_it); 
	     j_it != Cijk.j_end(k_it); ++j_it) {
	  OrdinalType j = index(j_it);
	  if (j >= k) {
	    ++j_coord_work[k_entry];
	    ++j_entry_count;
	  }
	}
      }
    }

    /*
    // Pad each row to have size divisible by alignment size
    enum { Align = KokkosArray::Impl::is_same<DeviceType,KokkosArray::Cuda>::value ? 32 : 2 };
    for ( size_type i = 0 ; i < dimension ; ++i ) {
      const size_t rem = coord_work[i] % Align;
      if (rem > 0) {
	const size_t pad = Align - rem;
	coord_work[i] += pad;
	entry_count += pad;
      }
    }
    */

    // Allocate tensor data
    FlatSparse3Tensor tensor ;
    tensor.m_k_coord = coord_array_type( "k_coord" , k_entry_count );
    tensor.m_j_coord = coord_array_type( "j_coord" , j_entry_count );
    tensor.m_value = value_array_type( "value" , j_entry_count );
    tensor.m_num_k = entry_array_type( "num_k" , dimension );
    tensor.m_num_j = entry_array_type( "num_j" , k_entry_count );
    tensor.m_k_row_map = row_map_array_type( "k_row_map" , dimension+1 );
    tensor.m_j_row_map = row_map_array_type( "j_row_map" , k_entry_count+1 );

    // Create mirror, is a view if is host memory
    typename coord_array_type::HostMirror
      host_k_coord = KokkosArray::create_mirror_view( tensor.m_k_coord );
    typename coord_array_type::HostMirror
      host_j_coord = KokkosArray::create_mirror_view( tensor.m_j_coord );
    typename value_array_type::HostMirror
      host_value = KokkosArray::create_mirror_view( tensor.m_value );
    typename entry_array_type::HostMirror
      host_num_k = KokkosArray::create_mirror_view( tensor.m_num_k );
    typename entry_array_type::HostMirror
      host_num_j = KokkosArray::create_mirror_view( tensor.m_num_j );
    typename entry_array_type::HostMirror
      host_k_row_map = KokkosArray::create_mirror_view( tensor.m_k_row_map );
    typename entry_array_type::HostMirror
      host_j_row_map = KokkosArray::create_mirror_view( tensor.m_j_row_map );

    // Compute k row map
    size_type sum = 0;
    host_k_row_map(0) = 0;
    for ( size_type i = 0 ; i < dimension ; ++i ) {
      sum += k_coord_work[i];
      host_k_row_map(i+1) = sum;
    }

    // Compute j row map
    sum = 0;
    host_j_row_map(0) = 0;
    for ( size_type i = 0 ; i < k_entry_count ; ++i ) {
      sum += j_coord_work[i];
      host_j_row_map(i+1) = sum;
    }

    for ( size_type i = 0 ; i < dimension ; ++i ) {
      k_coord_work[i] = host_k_row_map[i];
    }
    for ( size_type i = 0 ; i < k_entry_count ; ++i ) {
      j_coord_work[i] = host_j_row_map[i];
    }

    for (typename Cijk_type::i_iterator i_it=Cijk.i_begin(); 
	 i_it!=Cijk.i_end(); ++i_it) {
      OrdinalType i = index(i_it);
      for (typename Cijk_type::ik_iterator k_it = Cijk.k_begin(i_it); 
	   k_it != Cijk.k_end(i_it); ++k_it) {
	OrdinalType k = index(k_it);
	const size_type kEntry = k_coord_work[i]; 
	++k_coord_work[i];
	host_k_coord(kEntry) = k ;
	++host_num_k(i);
	for (typename Cijk_type::ikj_iterator j_it = Cijk.j_begin(k_it); 
	     j_it != Cijk.j_end(k_it); ++j_it) {
	  OrdinalType j = index(j_it);
	  ValueType c = Stokhos::value(j_it);
	  if (j >= k) {
	    const size_type jEntry = j_coord_work[kEntry]; 
	    ++j_coord_work[kEntry];
	    host_value(jEntry) = (j != k) ? c : 0.5*c;
	    host_j_coord(jEntry) = j ;
	    ++host_num_j(kEntry);
	    ++tensor.m_nnz;
	  }
	}
      }
    }

    // Copy data to device if necessary
    KokkosArray::deep_copy( tensor.m_k_coord , host_k_coord );
    KokkosArray::deep_copy( tensor.m_j_coord , host_j_coord );
    KokkosArray::deep_copy( tensor.m_value , host_value );
    KokkosArray::deep_copy( tensor.m_num_k , host_num_k );
    KokkosArray::deep_copy( tensor.m_num_j , host_num_j );
    KokkosArray::deep_copy( tensor.m_k_row_map , host_k_row_map );
    KokkosArray::deep_copy( tensor.m_j_row_map , host_j_row_map );

    tensor.m_flops = 5*tensor.m_nnz + dimension;

    return tensor ;
  }
};

template< class Device , typename OrdinalType , typename ValueType >
FlatSparse3Tensor<ValueType, Device>
create_flat_sparse_3_tensor( 
  const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk )
{
  return FlatSparse3Tensor<ValueType, Device>::create( basis, Cijk );
}

} /* namespace Stokhos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef STOKHOS_FLAT_SPARSE_3_TENSOR_HPP */



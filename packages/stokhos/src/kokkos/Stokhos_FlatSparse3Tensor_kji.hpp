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

#ifndef STOKHOS_FLAT_SPARSE_3_TENSOR_KJI_HPP
#define STOKHOS_FLAT_SPARSE_3_TENSOR_KJI_HPP

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
class FlatSparse3Tensor_kji {
public:

  typedef DeviceType                       device_type ;
  typedef typename device_type::size_type  size_type ;
  typedef ValueType                        value_type ;

private:

  typedef KokkosArray::View< size_type[] , device_type > coord_array_type ;
  typedef KokkosArray::View< value_type[], device_type > value_array_type ;
  typedef KokkosArray::View< size_type[], device_type > entry_array_type ;
  typedef KokkosArray::View< size_type[], device_type > row_map_array_type ;

  coord_array_type   m_j_coord ;
  coord_array_type   m_i_coord ;
  value_array_type   m_value ;
  entry_array_type   m_num_j ;
  entry_array_type   m_num_i ;
  row_map_array_type m_j_row_map ;
  row_map_array_type m_i_row_map ;
  size_type          m_nnz ;
  size_type          m_dim ;
  size_type          m_flops ;
 

public:

  inline
  ~FlatSparse3Tensor_kji() {}

  inline
  FlatSparse3Tensor_kji() : 
    m_j_coord() , 
    m_i_coord() , 
    m_value() , 
    m_num_j() , 
    m_num_i() , 
    m_j_row_map() , 
    m_i_row_map() , 
    m_nnz(0) ,
    m_dim(0) ,
    m_flops(0) {}

  inline
  FlatSparse3Tensor_kji( const FlatSparse3Tensor_kji & rhs ) : 
    m_j_coord( rhs.m_j_coord ) , 
    m_i_coord( rhs.m_i_coord ) ,
    m_value( rhs.m_value ) , 
    m_num_j( rhs.m_num_j ) ,
    m_num_i( rhs.m_num_i ) ,
    m_j_row_map( rhs.m_j_row_map ) , 
    m_i_row_map( rhs.m_i_row_map ) , 
    m_nnz( rhs.m_nnz ) ,
    m_dim( rhs.m_dim ) ,
    m_flops( rhs.m_flops ) {}

  inline
  FlatSparse3Tensor_kji & operator = ( const FlatSparse3Tensor_kji & rhs )
  {
    m_j_coord = rhs.m_j_coord ;
    m_i_coord = rhs.m_i_coord ;
    m_value = rhs.m_value ;
    m_num_j = rhs.m_num_j ;
    m_num_i = rhs.m_num_i ;
    m_j_row_map = rhs.m_j_row_map ;
    m_i_row_map = rhs.m_i_row_map ;
    m_nnz = rhs.m_nnz;
    m_dim = rhs.m_dim ;
    m_flops = rhs.m_flops ;
    return *this ;
  }

  /** \brief  Dimension of the tensor. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension() const { return m_dim ; }

  /** \brief  Number of k entries. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_k() const { return m_j_row_map.dimension(0) - 1 ; }

  /** \brief  Number of sparse entries. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type entry_count() const
  { return m_i_coord.dimension(0); }

  /** \brief  Begin j entries with a coordinate 'k' */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type j_begin( size_type k ) const
  { return m_j_row_map[k]; }

  /** \brief  End j entries with a coordinate 'k' */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type j_end( size_type k ) const
  { return m_j_row_map[k] + m_num_j(k); }

  /** \brief  Number of j entries with a coordinate 'k' */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_j( size_type k ) const
  { return m_num_j(k); }

  /** \brief  j coordinate for j entry 'jEntry' */
  KOKKOSARRAY_INLINE_FUNCTION
  const size_type& j_coord( const size_type jEntry ) const
  { return m_j_coord( jEntry ); }

  /** \brief  Begin i entries with a j entry 'jEntry' */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type i_begin( size_type jEntry ) const
  { return m_i_row_map[jEntry]; }

  /** \brief  End i entries with a j entry 'jEntry' */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type i_end( size_type jEntry ) const
  { return m_i_row_map[jEntry] + m_num_i(jEntry); }

  /** \brief  Number of i entries with a j entry 'jEntry' */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_i( size_type jEntry ) const
  { return m_num_i(jEntry); }

  /** \brief  i coordinate for i entry 'iEntry' */
  KOKKOSARRAY_INLINE_FUNCTION
  const size_type& i_coord( const size_type iEntry ) const
  { return m_i_coord( iEntry ); }

  /** \brief  Value for i entry 'iEntry' */
  KOKKOSARRAY_INLINE_FUNCTION
  const value_type & value( const size_type iEntry ) const
  { return m_value( iEntry ); }

  /** \brief Number of non-zero's */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_non_zeros() const 
  { return m_nnz; }

  /** \brief Number flop's per multiply-add */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_flops() const 
  { return m_flops; }

  template <typename OrdinalType>
  static FlatSparse3Tensor_kji
  create( const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
	  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk )
  {
    typedef Stokhos::Sparse3Tensor<OrdinalType,ValueType> Cijk_type;
    
    // Compute number of j's for each k
    const size_type dimension = basis.size();
    const size_type nk = Cijk.num_k();
    std::vector< size_t > j_coord_work( nk , (size_t) 0 );
    size_type j_entry_count = 0 ;
    for (typename Cijk_type::k_iterator k_it=Cijk.k_begin(); 
	 k_it!=Cijk.k_end(); ++k_it) {
      OrdinalType k = index(k_it);
      for (typename Cijk_type::kj_iterator j_it = Cijk.j_begin(k_it); 
	   j_it != Cijk.j_end(k_it); ++j_it) {
	OrdinalType j = index(j_it);
	if (j >= k) {
	  ++j_coord_work[k];
	  ++j_entry_count;
	}
      }
    }

    // Compute number of i's for each k and j
    std::vector< size_t > i_coord_work( j_entry_count , (size_t) 0 );
    size_type i_entry_count = 0 ;
    size_type j_entry = 0 ;
    for (typename Cijk_type::k_iterator k_it=Cijk.k_begin(); 
	 k_it!=Cijk.k_end(); ++k_it) {
      OrdinalType k = index(k_it);
      for (typename Cijk_type::kj_iterator j_it = Cijk.j_begin(k_it); 
	   j_it != Cijk.j_end(k_it); ++j_it) {
	OrdinalType j = index(j_it);
	if (j >= k) {
	  for (typename Cijk_type::kji_iterator i_it = Cijk.i_begin(j_it); 
	       i_it != Cijk.i_end(j_it); ++i_it) {
	    ++i_coord_work[j_entry];
	    ++i_entry_count;
	  }
	  ++j_entry;
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
    FlatSparse3Tensor_kji tensor ;
    tensor.m_dim = dimension;
    tensor.m_j_coord = coord_array_type( "j_coord" , j_entry_count );
    tensor.m_i_coord = coord_array_type( "i_coord" , i_entry_count );
    tensor.m_value = value_array_type( "value" , i_entry_count );
    tensor.m_num_j = entry_array_type( "num_j" , nk );
    tensor.m_num_i = entry_array_type( "num_i" , j_entry_count );
    tensor.m_j_row_map = row_map_array_type( "j_row_map" , nk+1 );
    tensor.m_i_row_map = row_map_array_type( "i_row_map" , j_entry_count+1 );
    tensor.m_flops = 3*j_entry_count + 2*i_entry_count;

    // Create mirror, is a view if is host memory
    typename coord_array_type::HostMirror
      host_j_coord = KokkosArray::create_mirror_view( tensor.m_j_coord );
    typename coord_array_type::HostMirror
      host_i_coord = KokkosArray::create_mirror_view( tensor.m_i_coord );
    typename value_array_type::HostMirror
      host_value = KokkosArray::create_mirror_view( tensor.m_value );
    typename entry_array_type::HostMirror
      host_num_j = KokkosArray::create_mirror_view( tensor.m_num_j );
    typename entry_array_type::HostMirror
      host_num_i = KokkosArray::create_mirror_view( tensor.m_num_i );
    typename entry_array_type::HostMirror
      host_j_row_map = KokkosArray::create_mirror_view( tensor.m_j_row_map );
    typename entry_array_type::HostMirror
      host_i_row_map = KokkosArray::create_mirror_view( tensor.m_i_row_map );

    // Compute j row map
    size_type sum = 0;
    host_j_row_map(0) = 0;
    for ( size_type k = 0 ; k < nk ; ++k ) {
      sum += j_coord_work[k];
      host_j_row_map(k+1) = sum;
    }

    // Compute i row map
    sum = 0;
    host_i_row_map(0) = 0;
    for ( size_type j = 0 ; j < j_entry_count ; ++j ) {
      sum += i_coord_work[j];
      host_i_row_map(j+1) = sum;
    }

    for ( size_type k = 0 ; k < nk ; ++k ) {
      j_coord_work[k] = host_j_row_map[k];
    }
    for ( size_type j = 0 ; j < j_entry_count ; ++j ) {
      i_coord_work[j] = host_i_row_map[j];
    }

    for (typename Cijk_type::k_iterator k_it=Cijk.k_begin(); 
	 k_it!=Cijk.k_end(); ++k_it) {
      OrdinalType k = index(k_it);
      for (typename Cijk_type::kj_iterator j_it = Cijk.j_begin(k_it); 
	   j_it != Cijk.j_end(k_it); ++j_it) {
	OrdinalType j = index(j_it);
	if (j >= k) {
	  const size_type jEntry = j_coord_work[k]; 
	  ++j_coord_work[k];
	  host_j_coord(jEntry) = j ;
	  ++host_num_j(k);
	  for (typename Cijk_type::kji_iterator i_it = Cijk.i_begin(j_it); 
	       i_it != Cijk.i_end(j_it); ++i_it) {
	    OrdinalType i = index(i_it);
	    ValueType c = Stokhos::value(i_it);
	    const size_type iEntry = i_coord_work[jEntry]; 
	    ++i_coord_work[jEntry];
	    host_value(iEntry) = (j != k) ? c : 0.5*c;
	    host_i_coord(iEntry) = i ;
	    ++host_num_i(jEntry);
	    ++tensor.m_nnz;
	  }
	}
      }
    }

    // Copy data to device if necessary
    KokkosArray::deep_copy( tensor.m_j_coord , host_j_coord );
    KokkosArray::deep_copy( tensor.m_i_coord , host_i_coord );
    KokkosArray::deep_copy( tensor.m_value , host_value );
    KokkosArray::deep_copy( tensor.m_num_j , host_num_j );
    KokkosArray::deep_copy( tensor.m_num_i , host_num_i );
    KokkosArray::deep_copy( tensor.m_j_row_map , host_j_row_map );
    KokkosArray::deep_copy( tensor.m_i_row_map , host_i_row_map );

    return tensor ;
  }
};

template< class Device , typename OrdinalType , typename ValueType >
FlatSparse3Tensor_kji<ValueType, Device>
create_flat_sparse_3_tensor_kji( 
  const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk )
{
  return FlatSparse3Tensor_kji<ValueType, Device>::create( basis, Cijk );
}

} /* namespace Stokhos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef STOKHOS_FLAT_SPARSE_3_TENSOR_KJI_HPP */



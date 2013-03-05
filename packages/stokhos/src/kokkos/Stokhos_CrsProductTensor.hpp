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

#ifndef STOKHOS_CRSPRODUCTTENSOR_HPP
#define STOKHOS_CRSPRODUCTTENSOR_HPP

#include "KokkosArray_View.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Stokhos {

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
template< typename ValueType , class DeviceType >
class CrsProductTensor {
public:

  typedef DeviceType                       device_type ;
  typedef typename device_type::size_type  size_type ;
  typedef ValueType                        value_type ;

private:

  typedef KokkosArray::View< value_type[] , device_type >  vec_type ;
  typedef KokkosArray::View< size_type[][2] , device_type > coord_array_type ;
  typedef KokkosArray::View< value_type[], device_type > value_array_type ;
  typedef KokkosArray::View< size_type[], device_type > entry_array_type ;
  typedef KokkosArray::View< size_type[], device_type > row_map_array_type ;

  coord_array_type   m_coord ;
  value_array_type   m_value ;
  entry_array_type   m_num_entry ;
  row_map_array_type m_row_map ;
  size_type          m_entry_max ;
  size_type          m_nnz ;
  size_type          m_flops ;

  struct CijkRowCount {
    unsigned count ;
    unsigned basis ;
    
    CijkRowCount()
      : count(0)
      , basis(0)
      {}
  };

public:

  inline
  ~CrsProductTensor() {}

  inline
  CrsProductTensor() : 
    m_coord() , m_value() , m_num_entry() , m_row_map() , 
    m_entry_max(0) , m_nnz(0) , m_flops(0) {}

  inline
  CrsProductTensor( const CrsProductTensor & rhs ) : 
    m_coord( rhs.m_coord ) , m_value( rhs.m_value ) , 
    m_num_entry( rhs.m_num_entry ) , m_row_map( rhs.m_row_map ) , 
    m_entry_max( rhs.m_entry_max ), m_nnz( rhs.m_nnz ), 
    m_flops( rhs.m_flops ) {}

  inline
  CrsProductTensor & operator = ( const CrsProductTensor & rhs )
  {
    m_coord = rhs.m_coord ;
    m_value = rhs.m_value ;
    m_num_entry = rhs.m_num_entry ;
    m_row_map = rhs.m_row_map ;
    m_entry_max = rhs.m_entry_max ;
    m_nnz = rhs.m_nnz;
    m_flops = rhs.m_flops;
    return *this ;
  }

  /** \brief  Dimension of the tensor. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension() const { return m_row_map.dimension(0) - 1 ; }

  /** \brief  Number of sparse entries. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type entry_count() const
  { return m_coord.dimension(0); }

  /** \brief  Maximum sparse entries for any coordinate */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type entry_maximum() const
  { return m_entry_max ; }

  /** \brief  Begin entries with a coordinate 'i' */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type entry_begin( size_type i ) const
  { return m_row_map[i]; }

  /** \brief  End entries with a coordinate 'i' */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type entry_end( size_type i ) const
  { return m_row_map[i] + m_num_entry(i); }

  /** \brief  Number of entries with a coordinate 'i' */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_entry( size_type i ) const
  { return m_num_entry(i); }

  /** \brief  Coordinates of an entry */
  KOKKOSARRAY_INLINE_FUNCTION
  const size_type& coord( const size_type entry , const size_type c ) const
  { return m_coord( entry , c ); }

  /** \brief  Value of an entry */
  KOKKOSARRAY_INLINE_FUNCTION
  const value_type & value( const size_type entry ) const
  { return m_value( entry ); }

  /** \brief Number of non-zero's */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_non_zeros() const 
  { return m_nnz; }

  /** \brief Number flop's per multiply-add */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_flops() const 
  { return m_flops; }

  template <typename OrdinalType>
  static CrsProductTensor
  create( const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
	  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk )
  {
    typedef Stokhos::Sparse3Tensor<OrdinalType,ValueType> Cijk_type;
    
    // Compute number of non-zeros for each i
    const size_type dimension = basis.size();
    std::vector< size_t > coord_work( dimension , (size_t) 0 );
    size_type entry_count = 0 ;
    for (typename Cijk_type::i_iterator i_it=Cijk.i_begin(); 
	 i_it!=Cijk.i_end(); ++i_it) {
      OrdinalType i = index(i_it);
      for (typename Cijk_type::ik_iterator k_it = Cijk.k_begin(i_it); 
	   k_it != Cijk.k_end(i_it); ++k_it) {
	OrdinalType k = index(k_it);
	for (typename Cijk_type::ikj_iterator j_it = Cijk.j_begin(k_it); 
	     j_it != Cijk.j_end(k_it); ++j_it) {
	  OrdinalType j = index(j_it);
	  if (j >= k) {
	    ++coord_work[i];
	    ++entry_count;
	  }
	}
      }
    }

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
    
    // Sort based on number of non-zeros
    std::vector< CijkRowCount > row_count( dimension );
    for ( size_type i = 0 ; i < dimension ; ++i ) {
      row_count[i].count = coord_work[i];
      row_count[i].basis = i;
    }
    //std::sort( row_count.begin() , row_count.end() , CompareCijkRowCount() );
    std::vector<size_type> sorted_row_map( dimension );
    for ( size_type i = 0 ; i < dimension ; ++i ) {
      coord_work[i] = row_count[i].count;
      sorted_row_map[ row_count[i].basis ] = i;
    }

    // Allocate tensor data
    CrsProductTensor tensor ;
    tensor.m_coord = coord_array_type( "tensor_coord" , entry_count );
    tensor.m_value = value_array_type( "tensor_value" , entry_count );
    tensor.m_num_entry = entry_array_type( "tensor_num_entry" , dimension );
    tensor.m_row_map = row_map_array_type( "tensor_row_map" , dimension+1 );
    tensor.m_entry_max = 0 ;

    // Create mirror, is a view if is host memory
    typename coord_array_type::HostMirror
      host_coord = KokkosArray::create_mirror_view( tensor.m_coord );
    typename value_array_type::HostMirror
      host_value = KokkosArray::create_mirror_view( tensor.m_value );
    typename entry_array_type::HostMirror
      host_num_entry = KokkosArray::create_mirror_view( tensor.m_num_entry );
    typename entry_array_type::HostMirror
      host_row_map = KokkosArray::create_mirror_view( tensor.m_row_map );

    // Compute row map
    size_type sum = 0;
    host_row_map(0) = 0;
    for ( size_type i = 0 ; i < dimension ; ++i ) {
      sum += coord_work[i];
      host_row_map(i+1) = sum;
    }

    for ( size_type iCoord = 0 ; iCoord < dimension ; ++iCoord ) {
      coord_work[iCoord] = host_row_map[iCoord];
    }

    for (typename Cijk_type::i_iterator i_it=Cijk.i_begin(); 
	 i_it!=Cijk.i_end(); ++i_it) {
      OrdinalType i = index(i_it);
      const size_type row = sorted_row_map[i];
      for (typename Cijk_type::ik_iterator k_it = Cijk.k_begin(i_it); 
	   k_it != Cijk.k_end(i_it); ++k_it) {
	OrdinalType k = index(k_it);
	for (typename Cijk_type::ikj_iterator j_it = Cijk.j_begin(k_it); 
	     j_it != Cijk.j_end(k_it); ++j_it) {
	  OrdinalType j = index(j_it);
	  ValueType c = Stokhos::value(j_it);
	  if (j >= k) {
	    const size_type n = coord_work[row]; ++coord_work[row];
	    host_value(n) = (j != k) ? c : 0.5*c;
	    host_coord(n,0) = j ;
	    host_coord(n,1) = k ;
	    ++host_num_entry(row);
	    ++tensor.m_nnz;
	  }
	}
      }
    }

    // Copy data to device if necessary
    KokkosArray::deep_copy( tensor.m_coord , host_coord );
    KokkosArray::deep_copy( tensor.m_value , host_value );
    KokkosArray::deep_copy( tensor.m_num_entry , host_num_entry );
    KokkosArray::deep_copy( tensor.m_row_map , host_row_map );

    for ( size_type i = 0 ; i < dimension ; ++i ) {
      tensor.m_entry_max = std::max( tensor.m_entry_max , host_num_entry(i) );
    }

    tensor.m_flops = 5*tensor.m_nnz + dimension;

    return tensor ;
  }
};

template< class Device , typename OrdinalType , typename ValueType >
CrsProductTensor<ValueType, Device>
create_product_tensor( 
  const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk )
{
  return CrsProductTensor<ValueType, Device>::create( basis, Cijk );
}

} /* namespace Stokhos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef STOKHOS_CRSPRODUCTTENSOR_HPP */



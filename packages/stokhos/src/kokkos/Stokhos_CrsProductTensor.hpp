// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_CRSPRODUCTTENSOR_HPP
#define STOKHOS_CRSPRODUCTTENSOR_HPP

#include "Kokkos_View.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Kokkos_Cuda.hpp"

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
template< typename ValueType, class DeviceType >
class CrsProductTensor {
public:

  typedef DeviceType                       device_type;
  typedef typename device_type::size_type  size_type;
  typedef ValueType                        value_type;

private:

  typedef Kokkos::View< value_type[], device_type >  vec_type;
  typedef Kokkos::View< size_type[], device_type > coord_array_type;
  typedef Kokkos::View< size_type[][2], device_type > coord2_array_type;
  typedef Kokkos::View< value_type[], device_type > value_array_type;
  typedef Kokkos::View< size_type[], device_type > entry_array_type;
  typedef Kokkos::View< size_type[], device_type > row_map_array_type;

  coord_array_type   m_coord;
  coord2_array_type  m_coord2;
  value_array_type   m_value;
  entry_array_type   m_num_entry;
  row_map_array_type m_row_map;
  size_type          m_entry_max;
  size_type          m_nnz;
  size_type          m_flops;

  struct CijkRowCount {
    unsigned count;
    unsigned basis;

    CijkRowCount()
      : count(0)
      , basis(0)
      {}
  };

  struct CompareCijkRowCount {
    bool operator() (const CijkRowCount& a, const CijkRowCount& b) const {
      return a.count < b.count;
    }
  };

public:

  inline
  ~CrsProductTensor() {}

  inline
  CrsProductTensor() :
    m_coord(),
    m_coord2(),
    m_value(),
    m_num_entry(),
    m_row_map(),
    m_entry_max(0),
    m_nnz(0),
    m_flops(0) {}

  inline
  CrsProductTensor( const CrsProductTensor & rhs ) :
    m_coord( rhs.m_coord ),
    m_coord2( rhs.m_coord2 ),
    m_value( rhs.m_value ),
    m_num_entry( rhs.m_num_entry ),
    m_row_map( rhs.m_row_map ),
    m_entry_max( rhs.m_entry_max ),
    m_nnz( rhs.m_nnz ),
    m_flops( rhs.m_flops ) {}

  inline
  CrsProductTensor & operator = ( const CrsProductTensor & rhs )
  {
    m_coord = rhs.m_coord;
    m_coord2 = rhs.m_coord2;
    m_value = rhs.m_value;
    m_num_entry = rhs.m_num_entry;
    m_row_map = rhs.m_row_map;
    m_entry_max = rhs.m_entry_max;
    m_nnz = rhs.m_nnz;
    m_flops = rhs.m_flops;
    return *this;
  }

  /** \brief  Dimension of the tensor. */
  KOKKOS_INLINE_FUNCTION
  size_type dimension() const { return m_row_map.dimension_0() - 1; }

  /** \brief  Number of sparse entries. */
  KOKKOS_INLINE_FUNCTION
  size_type entry_count() const
  { return m_coord.dimension_0(); }

  /** \brief  Maximum sparse entries for any coordinate */
  KOKKOS_INLINE_FUNCTION
  size_type entry_maximum() const
  { return m_entry_max; }

  /** \brief  Begin entries with a coordinate 'i' */
  KOKKOS_INLINE_FUNCTION
  size_type entry_begin( size_type i ) const
  { return m_row_map[i]; }

  /** \brief  End entries with a coordinate 'i' */
  KOKKOS_INLINE_FUNCTION
  size_type entry_end( size_type i ) const
  { return m_row_map[i] + m_num_entry(i); }

  /** \brief  Number of entries with a coordinate 'i' */
  KOKKOS_INLINE_FUNCTION
  size_type num_entry( size_type i ) const
  { return m_num_entry(i); }

  /** \brief  Coordinates of an entry */
  KOKKOS_INLINE_FUNCTION
  const size_type& coord( const size_type entry, const size_type c ) const
  { return m_coord2( entry, c ); }

  /** \brief  Coordinates of an entry */
  KOKKOS_INLINE_FUNCTION
  const size_type& coord( const size_type entry ) const
  { return m_coord( entry ); }

  /** \brief  Value of an entry */
  KOKKOS_INLINE_FUNCTION
  const value_type & value( const size_type entry ) const
  { return m_value( entry ); }

  /** \brief Number of non-zero's */
  KOKKOS_INLINE_FUNCTION
  size_type num_non_zeros() const
  { return m_nnz; }

  /** \brief Number flop's per multiply-add */
  KOKKOS_INLINE_FUNCTION
  size_type num_flops() const
  { return m_flops; }

  template <typename OrdinalType>
  static CrsProductTensor
  create( const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
          const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk,
          const Teuchos::ParameterList& params = Teuchos::ParameterList())
  {
    typedef Stokhos::Sparse3Tensor<OrdinalType,ValueType> Cijk_type;
    bool sort_nnz = false;
    if (params.isParameter("Sort Nonzeros"))
      sort_nnz = params.get<bool>("Sort Nonzeros");

    // Compute number of non-zeros for each i
    const size_type dimension = basis.size();
    std::vector< size_t > coord_work( dimension, (size_t) 0 );
    size_type entry_count = 0;
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
    enum { Align = Kokkos::Impl::is_same<DeviceType,Kokkos::Cuda>::value ? 32 : 2 };
    for ( size_type i = 0; i < dimension; ++i ) {
      const size_t rem = coord_work[i] % Align;
      if (rem > 0) {
        const size_t pad = Align - rem;
        coord_work[i] += pad;
        entry_count += pad;
      }
    }

    // Sort based on number of non-zeros
    std::vector< CijkRowCount > row_count( dimension );
    for ( size_type i = 0; i < dimension; ++i ) {
      row_count[i].count = coord_work[i];
      row_count[i].basis = i;
    }
    if (sort_nnz)
      std::sort( row_count.begin(), row_count.end(), CompareCijkRowCount() );
    std::vector<size_type> sorted_row_map( dimension );
    for ( size_type i = 0; i < dimension; ++i ) {
      coord_work[i] = row_count[i].count;
      sorted_row_map[ row_count[i].basis ] = i;
    }

    // Allocate tensor data
    CrsProductTensor tensor;
    tensor.m_coord = coord_array_type( "tensor_coord", entry_count );
    tensor.m_coord2 = coord2_array_type( "tensor_coord2", entry_count );
    tensor.m_value = value_array_type( "tensor_value", entry_count );
    tensor.m_num_entry = entry_array_type( "tensor_num_entry", dimension );
    tensor.m_row_map = row_map_array_type( "tensor_row_map", dimension+1 );
    tensor.m_entry_max = 0;

    // Create mirror, is a view if is host memory
    typename coord_array_type::HostMirror
      host_coord = Kokkos::create_mirror_view( tensor.m_coord );
    typename coord2_array_type::HostMirror
      host_coord2 = Kokkos::create_mirror_view( tensor.m_coord2 );
    typename value_array_type::HostMirror
      host_value = Kokkos::create_mirror_view( tensor.m_value );
    typename entry_array_type::HostMirror
      host_num_entry = Kokkos::create_mirror_view( tensor.m_num_entry );
    typename entry_array_type::HostMirror
      host_row_map = Kokkos::create_mirror_view( tensor.m_row_map );

    // Compute row map
    size_type sum = 0;
    host_row_map(0) = 0;
    for ( size_type i = 0; i < dimension; ++i ) {
      sum += coord_work[i];
      host_row_map(i+1) = sum;
      host_num_entry(i) = 0;
    }

    for ( size_type iCoord = 0; iCoord < dimension; ++iCoord ) {
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
            host_coord2(n,0) = j;
            host_coord2(n,1) = k;
            host_coord(n) = ( k << 16 ) | j;
            ++host_num_entry(row);
            ++tensor.m_nnz;
          }
        }
      }
    }

    // Copy data to device if necessary
    Kokkos::deep_copy( tensor.m_coord, host_coord );
    Kokkos::deep_copy( tensor.m_coord2, host_coord2 );
    Kokkos::deep_copy( tensor.m_value, host_value );
    Kokkos::deep_copy( tensor.m_num_entry, host_num_entry );
    Kokkos::deep_copy( tensor.m_row_map, host_row_map );

    for ( size_type i = 0; i < dimension; ++i ) {
      tensor.m_entry_max = std::max( tensor.m_entry_max, host_num_entry(i) );
    }

    tensor.m_flops = 5*tensor.m_nnz + dimension;

    return tensor;
  }
};

template< class Device, typename OrdinalType, typename ValueType >
CrsProductTensor<ValueType, Device>
create_product_tensor(
  const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk,
  const Teuchos::ParameterList& params = Teuchos::ParameterList())
{
  return CrsProductTensor<ValueType, Device>::create( basis, Cijk, params );
}

} /* namespace Stokhos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef STOKHOS_CRSPRODUCTTENSOR_HPP */

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

#ifndef STOKHOS_TILED_CRS_PRODUCT_TENSOR_HPP
#define STOKHOS_TILED_CRS_PRODUCT_TENSOR_HPP

#include "KokkosArray_View.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Stokhos_Sparse3TensorPartition.hpp"
#include "Teuchos_ParameterList.hpp"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Stokhos {

template< typename ValueType, class DeviceType >
class TiledCrsProductTensor {
public:

  typedef DeviceType                       device_type;
  typedef typename device_type::size_type  size_type;
  typedef ValueType                        value_type;

private:

  typedef KokkosArray::LayoutRight layout_type;
  typedef KokkosArray::View< value_type[], device_type >  vec_type;
  typedef KokkosArray::View< size_type[][2], device_type > coord_array_type;
  typedef KokkosArray::View< size_type[][3], device_type > coord_offset_type;
  typedef KokkosArray::View< size_type[][3], device_type > coord_range_type;
  typedef KokkosArray::View< value_type[], device_type > value_array_type;
  typedef KokkosArray::View< size_type**, layout_type, device_type > entry_array_type;
  typedef KokkosArray::View< size_type**, layout_type, device_type > row_map_array_type;
  typedef KokkosArray::View< size_type[], device_type > num_row_array_type;

  coord_array_type   m_coord;
  coord_offset_type  m_coord_offset;
  coord_range_type   m_coord_range;
  value_array_type   m_value;
  entry_array_type   m_num_entry;
  row_map_array_type m_row_map;
  num_row_array_type m_num_rows;
  size_type          m_dimension;
  size_type          m_tile_size;
  size_type          m_entry_max;
  size_type          m_max_num_rows;
  size_type          m_nnz;
  size_type          m_flops;

public:

  inline
  ~TiledCrsProductTensor() {}

  inline
  TiledCrsProductTensor() :
    m_coord(),
    m_coord_offset(),
    m_coord_range(),
    m_value(),
    m_num_entry(),
    m_row_map(),
    m_num_rows(),
    m_dimension(0),
    m_tile_size(0),
    m_entry_max(0),
    m_max_num_rows(0),
    m_nnz(0),
    m_flops(0) {}

  inline
  TiledCrsProductTensor( const TiledCrsProductTensor & rhs ) :
    m_coord( rhs.m_coord ),
    m_coord_offset( rhs.m_coord_offset ),
    m_coord_range( rhs.m_coord_range ),
    m_value( rhs.m_value ),
    m_num_entry( rhs.m_num_entry ),
    m_row_map( rhs.m_row_map ),
    m_num_rows( rhs.m_num_rows ),
    m_dimension( rhs.m_dimension ),
    m_tile_size( rhs.m_tile_size ),
    m_entry_max( rhs.m_entry_max ),
    m_max_num_rows( rhs.m_max_num_rows ),
    m_nnz( rhs.m_nnz ),
    m_flops( rhs.m_flops ) {}

  inline
  TiledCrsProductTensor & operator = ( const TiledCrsProductTensor & rhs )
  {
    m_coord = rhs.m_coord;
    m_coord_offset = rhs.m_coord_offset;
    m_coord_range = rhs.m_coord_range;
    m_value = rhs.m_value;
    m_num_entry = rhs.m_num_entry;
    m_row_map = rhs.m_row_map;
    m_num_rows = rhs.m_num_rows;
    m_dimension = rhs.m_dimension;
    m_tile_size = rhs.m_tile_size;
    m_entry_max = rhs.m_entry_max;
    m_max_num_rows = rhs.m_max_num_rows;
    m_nnz = rhs.m_nnz;
    m_flops = rhs.m_flops;
    return *this;
  }

  /** \brief  Dimension of the tensor. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension() const { return m_dimension; }

  /** \brief  Number of sparse entries. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type entry_count() const
  { return m_coord.dimension_0(); }

  /** \brief  Maximum sparse entries for any coordinate */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type entry_maximum() const
  { return m_entry_max; }

  /** \brief  Maximum number of rows in any tile. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type max_num_rows() const
  { return m_max_num_rows; }

  /** \brief  Number of rows in given tile. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_rows( size_type tile ) const
  { return m_num_rows(tile); }

  /** \brief  Begin entries with a coordinate 'i' */
  KOKKOSARRAY_INLINE_FUNCTION
  const size_type& entry_begin( size_type tile, size_type i ) const
  { return m_row_map(tile,i); }

  /** \brief  End entries with a coordinate 'i' */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type entry_end( size_type tile, size_type i ) const
  { return m_row_map(tile,i) + m_num_entry(tile,i); }

   /** \brief  Return row_map ptr */
  KOKKOSARRAY_INLINE_FUNCTION
  const size_type* row_map_ptr() const
  { return m_row_map.ptr_on_device(); }

  /** \brief  Number of entries with a coordinate 'i' */
  KOKKOSARRAY_INLINE_FUNCTION
  const size_type& num_entry( size_type tile, size_type i ) const
  { return m_num_entry(tile,i); }

  /** \brief  Coordinates of an entry */
  KOKKOSARRAY_INLINE_FUNCTION
  const size_type& coord( const size_type entry, const size_type c ) const
  { return m_coord( entry, c ); }

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

  /** \brief Number tiles */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type tile_size() const
  { return m_tile_size; }

  /** \brief Number tiles */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type num_tiles() const
  { return m_coord_offset.dimension_0(); }

  /** \brief Coordinate offset */
  KOKKOSARRAY_INLINE_FUNCTION
  const size_type& offset( const size_type entry, const size_type c ) const
  { return m_coord_offset( entry, c ); }

  /** \brief Coordinate range */
  KOKKOSARRAY_INLINE_FUNCTION
  const size_type& range( const size_type entry, const size_type c ) const
  { return m_coord_range( entry, c ); }

  template <typename OrdinalType>
  static TiledCrsProductTensor
  create( const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
          const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk,
          const Teuchos::ParameterList& params)
  {
    typedef Stokhos::Sparse3Tensor<OrdinalType,ValueType> Cijk_type;
    typedef Stokhos::CijkData<OrdinalType,ValueType> Cijk_Data_type;

    const size_type tile_size = params.get<int>("Tile Size");
    const size_type max_tiles = params.get<int>("Max Tiles");

    // Build tensor data list
    Teuchos::ArrayRCP<Cijk_Data_type> coordinate_list =
      Stokhos::build_cijk_coordinate_list(Cijk, Stokhos::CIJK_TWO_WAY_SYMMETRY);

    // Partition via RCB
    typedef Stokhos::RCB<Cijk_Data_type> rcb_type;
    typedef typename rcb_type::Box box_type;
    rcb_type rcb(tile_size, max_tiles, coordinate_list());
    Teuchos::RCP< Teuchos::Array< Teuchos::RCP<box_type> > > parts =
      rcb.get_parts();
    size_type num_parts = rcb.get_num_parts();

    // Compute number of non-zeros for each row in each part
    size_type total_num_rows = 0, max_num_rows = 0, entry_count = 0;
    Teuchos::Array< Teuchos::Array<size_type> > coord_work( num_parts );
    for (size_type part=0; part<num_parts; ++part) {
      Teuchos::RCP<box_type> box = (*parts)[part];
      size_type num_rows = box->delta_x;
      total_num_rows += num_rows;
      max_num_rows = std::max(max_num_rows, num_rows);
      coord_work[part].resize(num_rows, 0);

      size_type nc = box->coords.size();
      for (size_type c=0; c<nc; ++c) {
        size_type i = box->coords[c](0) - box->xmin;
        ++(coord_work[part][i]);
        ++entry_count;
      }
    }

    /*
    // Pad each row to have size divisible by alignment size
    enum { Align = KokkosArray::Impl::is_same<DeviceType,KokkosArray::Cuda>::value ? 32 : 2 };
    for ( size_type i = 0; i < dimension; ++i ) {
      const size_t rem = coord_work[i] % Align;
      if (rem > 0) {
        const size_t pad = Align - rem;
        coord_work[i] += pad;
        entry_count += pad;
      }
    }
    */

    // Allocate tensor data
    TiledCrsProductTensor tensor;
    tensor.m_coord =
      coord_array_type( "tensor_coord", entry_count );
    tensor.m_coord_offset =
      coord_offset_type( "tensor_coord_offset", num_parts );
    tensor.m_coord_range =
      coord_range_type( "tensor_coord_range", num_parts );
    tensor.m_value =
      value_array_type( "tensor_value", entry_count );
    tensor.m_num_entry =
      entry_array_type( "tensor_num_entry", num_parts, max_num_rows );
    tensor.m_row_map =
      row_map_array_type( "tensor_row_map", num_parts, max_num_rows+1 );
    tensor.m_num_rows =
      num_row_array_type( "tensor_num_rows", num_parts );
    tensor.m_dimension = basis.size();
    tensor.m_tile_size = tile_size;
    tensor.m_max_num_rows = max_num_rows;

    // Create mirror, is a view if is host memory
    typename coord_array_type::HostMirror host_coord =
      KokkosArray::create_mirror_view( tensor.m_coord );
    typename coord_offset_type::HostMirror host_coord_offset =
      KokkosArray::create_mirror_view( tensor.m_coord_offset );
    typename coord_range_type::HostMirror host_coord_range =
      KokkosArray::create_mirror_view( tensor.m_coord_range );
    typename value_array_type::HostMirror host_value =
      KokkosArray::create_mirror_view( tensor.m_value );
    typename entry_array_type::HostMirror host_num_entry =
      KokkosArray::create_mirror_view( tensor.m_num_entry );
    typename row_map_array_type::HostMirror host_row_map =
      KokkosArray::create_mirror_view( tensor.m_row_map );
    typename num_row_array_type::HostMirror host_num_rows =
      KokkosArray::create_mirror_view( tensor.m_num_rows );

    // Compute row map
    size_type sum = 0;
    for (size_type part=0; part<num_parts; ++part) {
      size_type nc = coord_work[part].size();
      host_row_map(part,0) = sum;
      for (size_type t=0; t<nc; ++t) {
        sum += coord_work[part][t];
        host_row_map(part,t+1) = sum;
      }
    }

    // Copy per part row offsets back into coord_work
    for (size_type part=0; part<num_parts; ++part) {
      size_type nc = coord_work[part].size();
      for (size_type t=0; t<nc; ++t) {
        coord_work[part][t] = host_row_map(part,t);
      }
    }

    // Fill in coordinate and value arrays
    for (size_type part=0; part<num_parts; ++part) {
      Teuchos::RCP<box_type> box = (*parts)[part];

      host_coord_offset(part,0) = box->xmin;
      host_coord_offset(part,1) = box->ymin;
      host_coord_offset(part,2) = box->zmin;

      host_coord_range(part,0) = box->delta_x;
      host_coord_range(part,1) = box->delta_y;
      host_coord_range(part,2) = box->delta_z;

      host_num_rows(part) = coord_work[part].size(); // also == box->delta_x

      size_type nc = box->coords.size();
      for (size_type t=0; t<nc; ++t) {
        const size_type i = box->coords[t].i;
        const size_type j = box->coords[t].j;
        const size_type k = box->coords[t].k;
        const value_type c = box->coords[t].c;

        const size_type row = i - box->xmin;
        const size_type n = coord_work[part][row];
        ++coord_work[part][row];

        host_value(n) = c;
        host_coord(n,0) = j - box->ymin;
        host_coord(n,1) = k - box->zmin;

        ++host_num_entry(part,row);
        ++tensor.m_nnz;
      }
    }

    // Copy data to device if necessary
    KokkosArray::deep_copy( tensor.m_coord, host_coord );
    KokkosArray::deep_copy( tensor.m_coord_offset, host_coord_offset );
    KokkosArray::deep_copy( tensor.m_coord_range, host_coord_range );
    KokkosArray::deep_copy( tensor.m_value, host_value );
    KokkosArray::deep_copy( tensor.m_num_entry, host_num_entry );
    KokkosArray::deep_copy( tensor.m_row_map, host_row_map );
    KokkosArray::deep_copy( tensor.m_num_rows, host_num_rows );

    tensor.m_entry_max = 0;
    tensor.m_flops = 0;
    for (size_type part=0; part<num_parts; ++part) {
      for ( size_type i = 0; i < host_num_rows(part); ++i ) {
        tensor.m_entry_max = std::max( tensor.m_entry_max,
                                       host_num_entry(part,i) );
        tensor.m_flops += 5*host_num_entry(part,i) + 1;
      }
    }

    return tensor;
  }
};

template< class Device, typename OrdinalType, typename ValueType >
TiledCrsProductTensor<ValueType, Device>
create_tiled_product_tensor(
  const Stokhos::ProductBasis<OrdinalType,ValueType>& basis,
  const Stokhos::Sparse3Tensor<OrdinalType,ValueType>& Cijk,
  const Teuchos::ParameterList& params)
{
  return TiledCrsProductTensor<ValueType, Device>::create(
    basis, Cijk, params );
}

} /* namespace Stokhos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef STOKHOS_TILED_CRS_PRODUCT_TENSOR_HPP */

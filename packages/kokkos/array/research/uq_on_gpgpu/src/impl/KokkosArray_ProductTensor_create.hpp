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

#ifndef KOKKOSARRAY_PRODUCTTENSOR_CREATE_HPP
#define KOKKOSARRAY_PRODUCTTENSOR_CREATE_HPP

#include <iostream>

#include <map>
#include <KokkosArray_CrsArray.hpp>

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  Create a sparse product tensor on the device
 *          from a map of tensor indices to values.
 *
 *  The std::map input guarantees uniqueness and proper sorting of
 *  the product tensor's symmetric entries.
 */
template< unsigned Rank , typename ValueType , class Device >
class CreateSparseProductTensor<
  SparseProductTensor< Rank , ValueType , Device > ,
  std::map< ProductTensorIndex<Rank> , ValueType > >
{
public:
  typedef SparseProductTensor<Rank,ValueType,Device> type ;
  typedef std::map< ProductTensorIndex< Rank > , ValueType > input_type ;

  static
  type create( const input_type & input )
  {
    typedef ValueType                   value_type ;
    typedef Device                      device_type ;
    typedef typename Device::size_type  size_type ;
    typedef typename type::map_type     map_type ;

    typedef View< value_type[], device_type> value_array_type ;

    type tensor ;

    tensor.m_coord = map_type( "sparse_tensor_coord" , input.size() );
    tensor.m_value = value_array_type( "sparse_tensor_value" , input.size() );

    // Try to create as a view, not a copy
    typename map_type::HostMirror
      host_coord = create_mirror_view( tensor.m_coord );

    // Try to create as a view, not a copy
    typename value_array_type::HostMirror
      host_value = create_mirror_view( tensor.m_value );

    size_type n = 0 ;

    tensor.m_dimen = 0 ;

    for ( typename input_type::const_iterator
          iter = input.begin() ; iter != input.end() ; ++iter , ++n ) {

      host_value(n) = (*iter).second ;

      for ( size_type c = 0 ; c < Rank ; ++c ) {

        const size_type i = (*iter).first.coord(c);

        host_coord(n,c) = i ;

        tensor.m_dimen = std::max( tensor.m_dimen , i + 1 );
      }
    }

    KokkosArray::deep_copy( tensor.m_coord , host_coord );
    KokkosArray::deep_copy( tensor.m_value , host_value );

    return tensor ;
  }
};

struct CijkRowCount {
  unsigned count ;
  unsigned basis ;

  CijkRowCount()
  : count(0)
  , basis(0)
  {}
};

struct CompareCijkRowCount {

  bool operator()( const CijkRowCount & lhs ,
                   const CijkRowCount & rhs ) const
  {
    return lhs.count != rhs.count ? lhs.count > rhs.count : (
      lhs.basis < rhs.basis );

  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  Create a sparse product tensor on the device
 *          from a map of tensor indices to values.
 *
 *  The std::map input guarantees uniqueness and proper sorting of
 *  the product tensor's symmetric entries.
 */
template< typename ValueType , class Device >
class CreateSparseProductTensor<
  CrsProductTensor< 3 , ValueType , Device > ,
  std::map< ProductTensorIndex<3> , ValueType > >
{
public:
  enum { Rank = 3 };
  typedef CrsProductTensor<Rank,ValueType,Device> type ;
  typedef std::map< ProductTensorIndex< Rank > , ValueType > input_type ;

  // input entries are sorted: coord(0) >= coord(1) >= coord(2)
  // thus last entry has maximum coordinate

  static
  type create( const input_type & input )
  {
    typedef ValueType                   value_type ;
    typedef Device                      device_type ;
    typedef typename Device::size_type  size_type ;

    //typedef CrsArray< size_type[2] , device_type > coord_array_type ;
    typedef View< size_type[][2] , device_type > coord_array_type ;
    typedef View< value_type[], device_type > value_array_type ;
    typedef View< size_type[], device_type > entry_array_type ;

    const size_type dimension =
      input.empty() ? 0 : 1 + (*input.rbegin()).first.coord(0);

    std::vector< size_t > coord_work( dimension , (size_t) 0 );

    size_type entry_count = 0 ;

    for ( typename input_type::const_iterator
          iter = input.begin() ; iter != input.end() ; ++iter ) {

      const size_type i = (*iter).first.coord(0);
      const size_type j = (*iter).first.coord(1);
      const size_type k = (*iter).first.coord(2);

      ++coord_work[i];
      ++entry_count ;
      if ( i != j ) { ++coord_work[j]; ++entry_count ; }
      if ( i != k && j != k ) { ++coord_work[k]; ++entry_count ; }
    }

    // Pad each row to have size divisble by 32
    enum { Align = Impl::is_same<Device,Cuda>::value ? 32 : 2 };
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

    type tensor ;

    //tensor.m_coord = create_crsarray< coord_array_type >( "tensor_coord" , coord_work );
    tensor.m_coord = coord_array_type( "tensor_coord" , entry_count );
    tensor.m_value = value_array_type( "tensor_value" , entry_count );
    tensor.m_num_entry = entry_array_type( "tensor_num_entry" , dimension );
    tensor.m_row_map = entry_array_type( "tensor_row_map" , dimension+1 );
    tensor.m_entry_max = 0 ;

/*
std::cout << std::endl << "CrsProductTensor" << std::endl
          << "  Tensor dimension     = " << dimension << std::endl
          << "  Tensor maximum entry = " << tensor.m_entry_max << std::endl
          << "  Compact tensor count = " << input.size() << std::endl
          << "  Crs     tensor count = " << entry_count << std::endl ;
*/

    // Create mirror, is a view if is host memory

    typename coord_array_type::HostMirror
      host_coord = create_mirror_view( tensor.m_coord );

    typename value_array_type::HostMirror
      host_value = create_mirror_view( tensor.m_value );

    typename entry_array_type::HostMirror
      host_num_entry = create_mirror_view( tensor.m_num_entry );

    typename entry_array_type::HostMirror
      host_row_map = create_mirror_view( tensor.m_row_map );

    // Fill arrays in coordinate order...

    size_type sum = 0;
    host_row_map(0) = 0;
    for ( size_type i = 0 ; i < dimension ; ++i ) {
      sum += coord_work[i];
      host_row_map(i+1) = sum;
    }

    for ( size_type iCoord = 0 ; iCoord < dimension ; ++iCoord ) {
      coord_work[iCoord] = host_row_map[iCoord];
    }

    for ( typename input_type::const_iterator
	    iter = input.begin() ; iter != input.end() ; ++iter ) {

      const size_type i = (*iter).first.coord(0);
      const size_type j = (*iter).first.coord(1);
      const size_type k = (*iter).first.coord(2);

      {
	const size_type row = sorted_row_map[i];
        const size_type n = coord_work[row]; ++coord_work[row];
        host_value(n) = (*iter).second ;
	if (j == k) host_value(n) *= 0.5;
        host_coord(n,0) = j ;
        host_coord(n,1) = k ;
	++host_num_entry(row);
	++tensor.m_nnz;
      }
      if ( i != j ) {
	const size_type row = sorted_row_map[j];
        const size_type n = coord_work[row]; ++coord_work[row];
        host_value(n) = (*iter).second ;
	if (i == k) host_value(n) *= 0.5;
        host_coord(n,0) = i ;
        host_coord(n,1) = k ;
	++host_num_entry(row);
	++tensor.m_nnz;
      }
      if ( i != k && j != k ) {
	const size_type row = sorted_row_map[k];
        const size_type n = coord_work[row]; ++coord_work[row];
        host_value(n) = (*iter).second ;
	if (i == j) host_value(n) *= 0.5;
        host_coord(n,0) = i ;
        host_coord(n,1) = j ;
	++host_num_entry(row);
	++tensor.m_nnz;
      }
    }

    // for (size_type i=0; i<dimension; ++i) {
    //   size_type iBeg = host_coord.row_map[i];
    //   size_type iEnd = host_coord.row_map[i+1];
    //   std::cout << "Row " << i << " has size " << iEnd-iBeg << std::endl;
    // }

    KokkosArray::deep_copy( tensor.m_coord , host_coord );
    KokkosArray::deep_copy( tensor.m_value , host_value );
    KokkosArray::deep_copy( tensor.m_num_entry , host_num_entry );
    KokkosArray::deep_copy( tensor.m_row_map , host_row_map );

    for ( size_type i = 0 ; i < dimension ; ++i ) {
      tensor.m_entry_max = std::max( tensor.m_entry_max , host_num_entry(i) );
    }

    return tensor ;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_PRODUCTTENSOR_CREATE_HPP */




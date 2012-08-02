/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_IMPL_CRSARRAY_FACTORY_HPP
#define KOKKOS_IMPL_CRSARRAY_FACTORY_HPP

#include <vector>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class > struct CrsArrayCreateMirror ;

template< class DataType , class LayoutType , class DeviceType , typename SizeType >
struct CrsArrayCreateMirror< CrsArray< DataType , LayoutType , DeviceType , SizeType > >
{
  typedef  CrsArray< DataType , LayoutType , DeviceType , SizeType > output_type ;

  inline static
  output_type create( const output_type & input ) { return input ; }

  template< class DeviceSrc >
  inline static
  output_type create( const CrsArray< DataType , LayoutType , DeviceSrc , SizeType > & input )
  {
    typedef View< SizeType[] , LayoutType , DeviceType > work_type ;

    work_type row_work( "mirror" , input.row_map.shape() );

    output_type output ;

    output.row_map = row_work ;
    output.entries = typename output_type::entries_type( "mirror" , input.entries.shape() );

    deep_copy( row_work ,       input.row_map );
    deep_copy( output.entries , input.entries );

    return output ;
  }
};

} // namespace Impl

template< class DataType , class LayoutType , class DeviceType , typename SizeType >
typename CrsArray< DataType , LayoutType , DeviceType , SizeType >::HostMirror
create_mirror_view( const CrsArray<DataType,LayoutType,DeviceType,SizeType > & input )
{
  typedef CrsArray< DataType , LayoutType , DeviceType , SizeType > input_type ;
  typedef typename input_type::HostMirror output_type ;

  return Impl::CrsArrayCreateMirror< output_type >::create( input );
}

template< class DataType , class LayoutType , class DeviceType , typename SizeType >
typename CrsArray< DataType , LayoutType , DeviceType , SizeType >::HostMirror
create_mirror( const CrsArray<DataType,LayoutType,DeviceType,SizeType > & input )
{
  typedef CrsArray< DataType , LayoutType , DeviceType , SizeType > input_type ;
  typedef typename input_type::HostMirror output_type ;

#if KOKKOS_MIRROR_VIEW_OPTIMIZE
  // Allow choice via type:
  return Impl::CrsArrayCreateMirror< output_type >::create( input );
#else
  // Force copy:
  return Impl::CrsArrayCreateMirror< output_type >::template create< DeviceType >( input );
#endif
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< class DataType ,
          class LayoutOutput ,
          class DeviceOutput ,
          typename MapSizeType ,
          typename InputSizeType >
struct Factory< CrsArray< DataType , LayoutOutput , DeviceOutput , MapSizeType > ,
                std::vector< InputSizeType > >
{
  typedef CrsArray< DataType , LayoutOutput , DeviceOutput, MapSizeType > output_type ;
  typedef std::vector< InputSizeType > input_type ;

  static
  output_type create( const std::string & label , const input_type & input )
  {
    typedef typename output_type::entries_type   entries_type ;
    typedef View< MapSizeType[] , LayoutOutput , DeviceOutput > work_type ;

    output_type output ;

    // Create the row map:

    const size_t length = input.size();

    {
      work_type row_work( "tmp" , length + 1 );

      typename work_type::HostMirror row_work_host =
        create_mirror_view( row_work );

      size_t sum = 0 ;
      row_work_host[0] = 0 ;
      for ( size_t i = 0 ; i < length ; ++i ) {
        row_work_host[i+1] = sum += input[i];
      }

      deep_copy( row_work , row_work_host );

      output.entries   = entries_type( label , sum );
      output.row_map   = row_work ;
    }

    return output ;
  }
};

//----------------------------------------------------------------------------

template< typename ValueType ,
          class LayoutOutput ,
          class DeviceOutput ,
          typename MapSizeType ,
          typename InputType >
struct Factory< CrsArray< ValueType , LayoutOutput , DeviceOutput , MapSizeType > ,
                std::vector< std::vector< InputType > > >
{
  typedef CrsArray< ValueType , LayoutOutput , DeviceOutput , MapSizeType > output_type ;
  typedef std::vector< std::vector< InputType > > input_type ;

  static const bool OK =
    StaticAssert< 1 == output_type::entries_type::Rank >::value ;

  static
  output_type create( const std::string & label , const input_type & input )
  {
    typedef typename output_type::entries_type   entries_type ;
    typedef View< MapSizeType[] , LayoutOutput , DeviceOutput > work_type ;

    output_type output ;

    // Create the row map:

    const size_t length = input.size();

    {
      work_type row_work( "tmp" , length + 1 );

      typename work_type::HostMirror row_work_host =
        create_mirror_view( row_work );

      size_t sum = 0 ;
      row_work_host[0] = 0 ;
      for ( size_t i = 0 ; i < length ; ++i ) {
        row_work_host[i+1] = sum += input[i].size();
      }

      deep_copy( row_work , row_work_host );

      output.entries   = entries_type( label , sum );
      output.row_map   = row_work ;
    }

    // Fill in the entries:
    {
      typename entries_type::HostMirror host_entries =
        create_mirror_view( output.entries );

      size_t sum = 0 ;
      for ( size_t i = 0 ; i < length ; ++i ) {
        for ( size_t j = 0 ; j < input[i].size() ; ++j , ++sum ) {
          host_entries( sum ) = input[i][j] ;
        }
      }

      deep_copy( output.entries , host_entries );
    }

    return output ;
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_CRSARRAY_FACTORY_HPP */


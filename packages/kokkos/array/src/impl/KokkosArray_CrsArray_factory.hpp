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
#include <impl/KokkosArray_MemoryView.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< class ArrayType , class Device , typename SizeType >
struct Factory< CrsArray< ArrayType , Device , SizeType > ,
                CrsArray< ArrayType , Device , SizeType > >
{
  typedef CrsArray< ArrayType, Device, SizeType > output_type ;

  static inline
  output_type create( const output_type & input )
  {
    typedef typename output_type::row_map_type row_map_type ;
    typedef typename output_type::entries_type entries_type ;

    output_type output ;

    output.row_map = Factory< row_map_type , row_map_type >
                       ::create( input.row_map );

    output.entries = Factory< entries_type , entries_type >
                       ::create( input.entries );

    Factory< entries_type , entries_type >
      ::deep_copy( output.entries , input.entries );

    return output ;
  }
};

//----------------------------------------------------------------------------

template< class ArrayType ,
          class DeviceDst , class DeviceSrc , typename SizeType >
struct Factory< CrsArray< ArrayType , DeviceDst , SizeType > ,
                CrsArray< ArrayType , DeviceSrc , SizeType > >
{
  typedef CrsArray< ArrayType, DeviceDst, SizeType > output_type ;
  typedef CrsArray< ArrayType, DeviceSrc, SizeType > input_type ;

  static inline
  output_type create( const input_type & input )
  {
    typedef typename output_type::row_map_type row_map_output_type ;
    typedef typename output_type::entries_type entries_output_type ;
    typedef typename input_type ::row_map_type row_map_input_type ;
    typedef typename input_type ::entries_type entries_input_type ;

    output_type output ;

    output.row_map = Factory< row_map_output_type , row_map_input_type >
                       ::create( input.row_map );

    output.entries = Factory< entries_output_type , entries_input_type >
                       ::create( input.entries );

    Factory< entries_output_type , entries_input_type >
      ::deep_copy( output.entries , input.entries );

    return output ;
  }
};

//----------------------------------------------------------------------------

template< class ArrayType , class DeviceOutput , typename SizeType >
struct Factory< CrsArray< ArrayType , DeviceOutput , SizeType > ,
                MirrorUseView >
{
  typedef CrsArray< ArrayType , DeviceOutput , SizeType > output_type ;

  static inline
  const output_type & create( const output_type & input ) { return input ; }

  template< class DeviceInput >
  static inline
  output_type create(
    const CrsArray< ArrayType , DeviceInput , SizeType > & input )
  {
    typedef typename output_type::row_map_type row_map_type ;
    typedef typename output_type::entries_type entries_type ;

    output_type output ;
    output.row_map = Factory< row_map_type , MirrorUseView >
                       ::create( input.row_map );
    output.entries = Factory< entries_type , MirrorUseView >
                       ::create( input.entries );

    return output ;
  }
};

//----------------------------------------------------------------------------

template< class ArrayType ,
          class DeviceOutput ,
          typename MapSizeType ,
          typename InputSizeType >
struct Factory< CrsArray< ArrayType , DeviceOutput , MapSizeType > ,
                std::vector< InputSizeType > >
{
  typedef CrsArray< ArrayType , DeviceOutput, MapSizeType > output_type ;
  typedef std::vector< InputSizeType > input_type ;

  static
  output_type create( const std::string & label , const input_type & input )
  {
    typedef typename output_type::row_map_type row_map_type ;
    typedef typename output_type::entries_type entries_type ;

    output_type output ;

    output.row_map = Factory< row_map_type , input_type >
                       ::create( label , input );

    output.entries = Factory< entries_type , void >
                       ::create( label , output.row_map.sum() );

    return output ;
  }
};

//----------------------------------------------------------------------------

template< typename ValueType ,
          class DeviceOutput ,
          typename MapSizeType ,
          typename InputType >
struct Factory< CrsArray< ValueType , DeviceOutput , MapSizeType > ,
                std::vector< std::vector< InputType > > >
{
  typedef CrsArray< ValueType , DeviceOutput , MapSizeType > output_type ;
  typedef std::vector< std::vector< InputType > > input_type ;

  static const bool OK =
    StaticAssert< 1 == output_type::entries_type::Rank >::value ;

  static
  output_type create( const std::string & label , const input_type & input )
  {
    typedef typename output_type::row_map_type   row_map_type ;
    typedef typename output_type::entries_type   entries_type ;
    typedef typename DeviceOutput::memory_space  memory_space ;

    typedef MemoryView< MapSizeType , memory_space > data_type ;
    typedef typename data_type::HostMirror data_mirror_type ;

    // Create the row map:

    const size_t length = input.size();

    output_type output ;

    output.row_map.m_length = length ;
    output.row_map.m_data.allocate( length + 1 , label );

    data_mirror_type tmp = Factory< data_mirror_type , MirrorUseView >
                             ::create( output.row_map.m_data , length + 1 );

    output.row_map.m_sum = tmp[0] = 0 ;
    for ( size_t i = 0 ; i < length ; ++i ) {
      tmp[i+1] = output.row_map.m_sum += input[i].size();
    }

    Factory< data_type , data_mirror_type >
      ::deep_copy( output.row_map.m_data , tmp , length + 1 );

    // Create an populate the entries:

    output.entries = Factory< entries_type , void >
                       ::create( label , output.row_map.m_sum );

    typename entries_type::HostMirror host_entries =
      Factory< typename entries_type::HostMirror , MirrorUseView >
        ::create( output.entries );

    size_t total_count = 0 ;
    for ( size_t i = 0 ; i < length ; ++i ) {
      for ( size_t j = 0 ; j < input[i].size() ; ++j , ++total_count ) {
        host_entries( total_count ) = input[i][j] ;
      }
    }

    Factory< entries_type , typename entries_type::HostMirror >
      ::deep_copy( output.entries , host_entries );

    return output ;
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_CRSARRAY_FACTORY_HPP */


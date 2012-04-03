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
#include <impl/Kokkos_MemoryView.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class ArrayType , class Device , typename SizeType >
struct Factory< CrsArray< ArrayType , Device , SizeType > ,
                CrsArray< ArrayType , Device , SizeType > >
{
  typedef CrsArray< ArrayType, Device, SizeType > output_type ;
  typedef SizeType                          size_type ;
  typedef typename output_type::value_type  value_type ;
  typedef typename Device::memory_space     device_memory ;

  typedef MemoryView< size_type , device_memory >  output_row ;
  typedef MemoryView< value_type, device_memory >  output_data ;

  static inline
  void deep_copy( const output_type & output , const output_type & input )
  {
    Factory< output_row , output_row >
      ::deep_copy( output.m_row_map, input.m_row_map , output.m_row_count + 1 );

    Factory< output_data , output_data >
      ::deep_copy( output.m_data, input.m_data , output.m_index_map.allocation_size() );
  }

  static inline
  output_type create( const output_type & input )
  {
    output_type output ;
    output.m_row_count = input.m_row_count ;
    output.m_index_map = input.m_index_map ;
    output.m_row_map.allocate( output.m_row_count + 1 , std::string() );
    output.m_data.allocate( output.m_index_map.allocation_size() , std::string() );

    deep_copy( output , input );

    return output ;
  }
};

template< class Device , typename SizeType >
struct Factory< CrsArray< void , Device , SizeType > ,
                CrsArray< void , Device , SizeType > >
{
  typedef CrsArray< void, Device, SizeType > output_type ;

  static inline
  output_type create( const output_type & input )
  {
    typedef SizeType                          size_type ;
    typedef typename Device::memory_space  device_memory ;

    typedef MemoryView< size_type , device_memory >  output_row ;

    output_type output ;
    output.m_row_count = input.m_row_count ;
    output.m_row_map.allocate( output.m_row_count + 1 , std::string() );

    Factory< output_row , output_row >
      ::deep_copy( output.m_row_map, input.m_row_map , output.m_row_count + 1 );

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
    typedef CrsArray< ArrayType , DeviceInput , SizeType > input_type ;
    return Factory< output_type , input_type >::create( input );
  }
};

//----------------------------------------------------------------------------

template< class DeviceOutput , typename SizeType >
struct Factory< CrsArray< void , DeviceOutput , SizeType > ,
                MirrorUseView >
{
  typedef CrsArray< void , DeviceOutput , SizeType > output_type ;

  static inline
  const output_type & create( const output_type & input ) { return input ; }

  template< class DeviceInput >
  static inline
  output_type create(
    const CrsArray< void , DeviceInput , SizeType > & input )
  {
    typedef CrsArray< void , DeviceInput , SizeType > input_type ;
    return Factory< output_type , input_type >::create( input );
  }
};

//----------------------------------------------------------------------------

template< class DeviceOutput , typename MapSizeType , typename InputSizeType >
struct Factory<
        CrsArray< void , DeviceOutput , MapSizeType > ,
        std::vector< InputSizeType > >
{
  typedef CrsArray< void , DeviceOutput, MapSizeType > output_type ;
  typedef std::vector< InputSizeType > input_type ;

  static
  output_type create( const std::string & label , const input_type & input )
  {
    typedef MemoryView< MapSizeType , typename DeviceOutput::memory_space > memory_output ;
    typedef typename memory_output::HostMirror  memory_mirror ;

    const size_t row_count = input.size();

    output_type output ;

    output.m_row_map.allocate( row_count + 1 , label );
    output.m_row_count = row_count ;

    // If same memory space then a view:
    memory_mirror tmp = Factory< memory_mirror , MirrorUseView >
                          ::create( output.m_row_map , row_count + 1 );

    tmp[0] = 0 ;
    for ( size_t i = 0 ; i < row_count ; ++i ) {
      tmp[i+1] = tmp[i] + input[i] ;
    }

    Factory< memory_output , memory_mirror >
      ::deep_copy( output.m_row_map , tmp , row_count + 1 );

    return output ;
  }
};

//----------------------------------------------------------------------------

template< class ArrayType ,
          class DeviceOutput ,
          typename MapSizeType ,
          typename InputSizeType >
struct Factory<
        CrsArray< ArrayType , DeviceOutput , MapSizeType > ,
        std::vector< InputSizeType > >
{
  typedef CrsArray< ArrayType , DeviceOutput, MapSizeType > output_type ;
  typedef std::vector< InputSizeType > input_type ;

  static
  output_type create( const std::string & label , const input_type & input )
  {
    typedef typename output_type::value_type          value_type ;
    typedef typename DeviceOutput::memory_space       output_space ;
    typedef MemoryView< MapSizeType , output_space >  output_row_type ;
    typedef typename output_row_type::HostMirror      mirror_row_type ;

    const size_t row_count = input.size();

    output_type output ;

    output.m_row_map.allocate( row_count + 1 , label );
    output.m_row_count = row_count ;

    // If same memory space then a view:
    mirror_row_type tmp =
      Factory< mirror_row_type , MirrorUseView >
        ::create( output.m_row_map , row_count + 1 );

    tmp[0] = 0 ;
    for ( size_t i = 0 ; i < row_count ; ++i ) {
      tmp[i+1] = tmp[i] + input[i] ;
    }

    output.m_index_map.template assign<value_type>( tmp[row_count] );
    output.m_data.allocate( output.m_index_map.allocation_size() , label );

    Factory< output_row_type , mirror_row_type >
      ::deep_copy( output.m_row_map , tmp , row_count + 1 );

    return output ;
  }
};

//----------------------------------------------------------------------------

template< class ArrayType ,
          class DeviceOutput ,
          typename MapSizeType ,
          typename InputType >
struct Factory< CrsArray< ArrayType , DeviceOutput , MapSizeType > ,
                std::vector< std::vector< InputType > > >
{
  typedef CrsArray< ArrayType , DeviceOutput , MapSizeType > output_type ;
  typedef std::vector< std::vector< InputType > > input_type ;

  static
  output_type create( const std::string & label , const input_type & input )
  {
    typedef typename DeviceOutput::memory_space memory_output ;
    typedef typename output_type::value_type    value_type ;

    typedef MemoryView< MapSizeType , memory_output > row_type ;
    typedef MemoryView< value_type ,  memory_output > data_type ;
    typedef MemoryView< MapSizeType , Host >          host_row_type ;
    typedef MemoryView< value_type ,  Host >          host_data_type ;

    const size_t row_count = input.size();
    size_t total_count = 0 ;
    for ( size_t i = 0 ; i < row_count ; ++i ) {
      total_count += input[i].size();
    }

    output_type crs ;

    crs.m_row_count = row_count ;
    crs.m_row_map.allocate( row_count + 1 , label );
    crs.m_index_map.template assign<value_type>( total_count );
    crs.m_data.allocate( crs.m_index_map.allocation_size() , label );

    // If same memory space then a view:
    host_row_type host_row =
      Factory< host_row_type , MirrorUseView >
        ::create( crs.m_row_map , row_count + 1 );

    host_data_type host_data =
      Factory< host_data_type , MirrorUseView >
        ::create( crs.m_data , crs.m_index_map.allocation_size() );

    host_row[0] = 0 ;
    for ( size_t i = 0 ; i < row_count ; ++i ) {
      host_row[i+1] = host_row[i] + input[i].size();
    }

    Factory< row_type , host_row_type >
      ::deep_copy( crs.m_row_map , host_row , row_count + 1 );

    switch( crs.entry_rank() ) {
    case 1 :
      total_count = 0 ;
      for ( size_t i = 0 ; i < row_count ; ++i ) {
        for ( size_t j = 0 ; j < input[i].size() ; ++j , ++total_count ) {
          host_data[ total_count ] = input[i][j] ;
        }
      }
      break ;
    }

    Factory< data_type , host_data_type >
      ::deep_copy( crs.m_data , host_data , total_count );

    return crs ;
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_CRSARRAY_FACTORY_HPP */


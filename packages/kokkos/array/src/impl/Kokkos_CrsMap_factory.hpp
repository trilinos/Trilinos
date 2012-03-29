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

#ifndef KOKKOS_IMPL_CRSMAP_FACTORY_HPP
#define KOKKOS_IMPL_CRSMAP_FACTORY_HPP

#include <vector>
#include <impl/Kokkos_MemoryView.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {



} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class DeviceOutput , class DeviceInput , typename SizeType >
struct Factory< CrsMap< DeviceOutput, CrsColumnMap , SizeType > ,
                CrsMap< DeviceInput , CrsColumnMap , SizeType > >
{
  typedef CrsMap< DeviceOutput, CrsColumnMap , SizeType > output_type ;
  typedef CrsMap< DeviceInput,  CrsColumnMap , SizeType > input_type ;

  static inline
  output_type create( const input_type & input )
  {
    typedef MemoryView< SizeType , typename DeviceOutput::memory_space > memory_output ;
    typedef MemoryView< SizeType , typename DeviceInput::memory_space >  memory_input ;

    const size_t offset_count = input.m_row_count + 1 + input.m_entry_count ;
    
    output_type output ;
    output.m_row_count    = input.m_row_count ;
    output.m_entry_count  = input.m_entry_count ;

    output.m_memory.allocate( offset_count , std::string() );
    output.m_column.m_map =
      output.m_memory.ptr_on_device() + output.m_row_count + 1 ;

    Factory< memory_output , memory_input >
      ::deep_copy( output.m_memory, input.m_memory , offset_count );

    return output ;
  }
};

//----------------------------------------------------------------------------

template< class DeviceOutput , typename SizeType >
struct Factory< CrsMap< DeviceOutput , CrsColumnMap , SizeType > ,
                MirrorUseView >
{
  typedef CrsMap< DeviceOutput, CrsColumnMap , SizeType > output_type ;

  static inline
  const output_type & create( const output_type & input ) { return input ; }

  template< class DeviceInput >
  static inline
  output_type create(
    const CrsMap< DeviceInput , CrsColumnMap , SizeType > & input )
  {
    typedef CrsMap< DeviceInput , CrsColumnMap , SizeType > input_type ;
    return Factory< output_type , input_type >::create( input );
  }
};

//----------------------------------------------------------------------------

template< class DeviceOutput , class DeviceInput , typename SizeType >
struct Factory< CrsMap< DeviceOutput, CrsColumnIdentity , SizeType > ,
                CrsMap< DeviceInput , CrsColumnIdentity , SizeType > >
{
  typedef CrsMap< DeviceOutput, CrsColumnIdentity , SizeType > output_type ;
  typedef CrsMap< DeviceInput , CrsColumnIdentity , SizeType > input_type ;

  static inline
  output_type create( const input_type & input )
  {
    typedef MemoryView< SizeType , typename DeviceOutput::memory_space > memory_output ;
    typedef MemoryView< SizeType , typename DeviceInput::memory_space >  memory_input ;

    output_type output ;
    output.m_row_count    = input.m_row_count ;
    output.m_entry_count  = input.m_entry_count ;

    output.m_memory.allocate( output.m_row_count + 1 , std::string() );

    Factory< memory_output , memory_input >
      ::deep_copy( output.m_memory, input.m_memory , output.m_row_count + 1 );

    return output ;
  }
};

//----------------------------------------------------------------------------

template< class DeviceOutput , typename SizeType >
struct Factory< CrsMap< DeviceOutput , CrsColumnIdentity , SizeType > ,
                MirrorUseView >
{
  typedef CrsMap< DeviceOutput, CrsColumnIdentity , SizeType > output_type ;

  static inline
  const output_type & create( const output_type & input ) { return input ; }

  template< class DeviceInput >
  static inline
  output_type create(
    const CrsMap< DeviceInput , CrsColumnIdentity , SizeType > & input )
  {
    typedef CrsMap< DeviceInput , CrsColumnIdentity , SizeType > input_type ;
    return Factory< output_type , input_type >::create( input );
  }
};

//----------------------------------------------------------------------------

template< class DeviceOutput , typename MapSizeType , typename InputSizeType >
struct Factory<
        CrsMap< DeviceOutput , CrsColumnIdentity , MapSizeType > ,
        std::vector< InputSizeType > >
{
  typedef CrsMap< DeviceOutput, CrsColumnIdentity, MapSizeType > output_type ;
  typedef std::vector< InputSizeType > input_type ;

  static
  output_type create( const std::string & label , const input_type & input )
  {
    typedef MemoryView< MapSizeType , typename DeviceOutput::memory_space > memory_output ;
    typedef typename memory_output::HostMirror  memory_mirror ;

    const size_t row_count = input.size();

    output_type output ;

    output.m_memory.allocate( row_count + 1 , label );
    output.m_row_count    = row_count ;
    output.m_entry_count  = 0 ;

    // If same memory space then a view:
    memory_mirror tmp = Factory< memory_mirror , MirrorUseView >
                          ::create( output.m_memory , row_count + 1 );

    tmp[0] = 0 ;
    for ( size_t i = 0 ; i < row_count ; ++i ) {
      tmp[i+1] = output.m_entry_count += input[i] ;
    }

    Factory< memory_output , memory_mirror >
      ::deep_copy( output.m_memory , tmp , row_count + 1 );

    return output ;
  }
};

//----------------------------------------------------------------------------

template< class DeviceOutput , typename MapSizeType , typename InputSizeType >
struct Factory< CrsMap< DeviceOutput , CrsColumnMap , MapSizeType > ,
                std::vector< std::vector< InputSizeType > > >
{
  typedef CrsMap< DeviceOutput , CrsColumnMap , MapSizeType > output_type ;
  typedef std::vector< std::vector< InputSizeType > > input_type ;

  static
  output_type create( const std::string & label , const input_type & input )
  {
    typedef MemoryView< MapSizeType , typename DeviceOutput::memory_space > memory_output ;
    typedef typename memory_output::HostMirror  memory_mirror ;

    const size_t row_count = input.size();
    size_t total_count = 0 ;
    for ( size_t i = 0 ; i < row_count ; ++i ) {
      total_count += input[i].size();
    }
    // Total number of members required:
    size_t offset_count = total_count + row_count + 1 ;

    output_type crs ;

    crs.m_memory.allocate( offset_count , label );
    crs.m_column.m_map = crs.m_memory.ptr_on_device() + row_count + 1 ;
    crs.m_row_count    = row_count ;
    crs.m_entry_count  = total_count ;

    // If same memory space then a view:
    memory_mirror tmp = Factory< memory_mirror , MirrorUseView >
                          ::create( crs.m_memory , offset_count );

    tmp[0] = 0 ;
    total_count = 0 ;
    offset_count = row_count + 1 ;
    for ( size_t i = 0 ; i < row_count ; ++i ) {
      const std::vector< InputSizeType > & column = input[i] ;
      for ( size_t j = 0 ; j < column.size() ; ++j , ++offset_count ) {
        tmp[offset_count] = column[j] ;
      }
      tmp[i+1] = total_count += column.size();
    }

    Factory< memory_output , memory_mirror >
      ::deep_copy( crs.m_memory , tmp , offset_count );

    return crs ;
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_CRSMAP_FACTORY_HPP */


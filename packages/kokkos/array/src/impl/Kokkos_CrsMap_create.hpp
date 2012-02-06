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

#ifndef KOKKOS_IMPL_CRSMAP_CREATE_HPP
#define KOKKOS_IMPL_CRSMAP_CREATE_HPP

#include <vector>
#include <impl/Kokkos_MemoryView.hpp>

namespace Kokkos {
namespace Impl {

template< class , class > class CreateCrsMap ;

//----------------------------------------------------------------------------

template< class Device ,
          template< class , typename > class ColumnType ,
          typename SizeType >
class CreateMirror< CrsMap< Device , ColumnType , SizeType > , true > {
public:
  typedef CrsMap< Device , ColumnType , SizeType >     view_type ;
  typedef typename view_type::HostMirror  mirror_type ;

  static
  mirror_type create( const view_type & v )
    { return v ; }
};

template< class DeviceSrc , typename SizeType >
class CreateCrsMap< CrsMap<Host,CrsColumnMap,SizeType> ,
                    CrsMap<DeviceSrc,CrsColumnMap,SizeType> > {
public:
  typedef CrsMap< Host , CrsColumnMap , SizeType >       type ;
  typedef CrsMap< DeviceSrc , CrsColumnMap , SizeType >  input_type ;

  static
  type create( const input_type & input )
  {
    typedef typename DeviceSrc::memory_space                memory_space ;
    typedef          MemoryView< SizeType , memory_space >  memory_type ;
    typedef typename memory_type::HostMirror                mirror_type ;

    const size_t offset_count = input.m_row_count + 1 + input.m_entry_count ;
    
    type crs ;
    crs.m_memory.allocate( offset_count , std::string() );
    crs.m_column.m_map = crs.m_memory.ptr_on_device() + input.m_row_count + 1 ;
    crs.m_row_count    = input.m_row_count ;
    crs.m_entry_count  = input.m_entry_count ;

    DeepCopy< mirror_type , memory_type >
      ::run( crs.m_memory, input.m_memory , offset_count );

    return crs ;
  }
};

template< class DeviceSrc , typename SizeType >
class CreateCrsMap< CrsMap<Host,CrsColumnIdentity,SizeType> ,
                    CrsMap<DeviceSrc,CrsColumnIdentity,SizeType> > {
public:
  typedef CrsMap< Host , CrsColumnIdentity , SizeType >       type ;
  typedef CrsMap< DeviceSrc , CrsColumnIdentity , SizeType >  input_type ;

  static
  type create( const input_type & input )
  {
    typedef typename DeviceSrc::memory_space                memory_space ;
    typedef          MemoryView< SizeType , memory_space >  memory_type ;
    typedef typename memory_type::HostMirror                mirror_type ;

    const size_t offset_count = input.m_row_count + 1 ;
    
    type crs ;
    crs.m_memory.allocate( offset_count , std::string() );
    crs.m_row_count   = input.m_row_count ;
    crs.m_entry_count = input.m_entry_count ;

    DeepCopy< mirror_type , memory_type >
      ::run( crs.m_memory, input.m_memory , offset_count );

    return crs ;
  }
};

template< class Device ,
          template< class , typename > class ColumnType ,
          typename SizeType >
class CreateMirror< CrsMap< Device , ColumnType , SizeType > , false > {
public:
  typedef CrsMap< Device , ColumnType , SizeType >     view_type ;
  typedef typename view_type::HostMirror  mirror_type ;

  static
  mirror_type create( const view_type & v )
    { return CreateCrsMap< mirror_type , view_type >::create( v ); }
};

//----------------------------------------------------------------------------

template< class Device , typename SizeType >
class CreateCrsMap< Device , std::vector< SizeType > > {
public:

  typedef CrsMap< Device , CrsColumnIdentity > type ;

  static
  type create( const std::string & label ,
               const std::vector< SizeType > & counts )
  {
    typedef typename Device::size_type                       size_type ;
    typedef typename Device::memory_space                    memory_space ;
    typedef          MemoryView< size_type , memory_space >  memory_type ;
    typedef typename memory_type::HostMirror                 mirror_type ;

    enum { is_host_memory = SameType< memory_space , Host >::value };

    typedef CreateMirror< memory_type , is_host_memory > create_mirror ;

    const size_t row_count = counts.size();

    type crs ;

    crs.m_memory.allocate( row_count + 1 , label );
    crs.m_row_count    = row_count ;
    crs.m_entry_count  = 0 ;

    // If same memory space then a view:
    mirror_type tmp = create_mirror::create( crs.m_memory , row_count + 1 );

    tmp[0] = 0 ;
    for ( size_t i = 0 ; i < row_count ; ++i ) {
      tmp[i+1] = crs.m_entry_count += counts[i] ;
    }

    DeepCopy< memory_type , mirror_type >::run( crs.m_memory , tmp , row_count + 1 );

    return crs ;
  }
};

//----------------------------------------------------------------------------

template< class Device , typename SizeType >
class CreateCrsMap< Device , std::vector< std::vector< SizeType > > > {
public:

  typedef CrsMap< Device , CrsColumnMap > type ;

  static
  type create( const std::string & label ,
               const std::vector< std::vector< SizeType > > & input )
  {
    typedef typename Device::size_type                       size_type ;
    typedef typename Device::memory_space                    memory_space ;
    typedef          MemoryView< size_type , memory_space >  memory_type ;
    typedef typename memory_type::HostMirror                 mirror_type ;

    enum { is_host_memory = SameType< memory_space , Host >::value };

    typedef CreateMirror< memory_type , is_host_memory > create_mirror ;

    const size_t row_count = input.size();
    size_type total_count = 0 ;
    for ( size_t i = 0 ; i < row_count ; ++i ) {
      total_count += input[i].size();
    }
    // Total number of members required:
    size_type offset_count = total_count + row_count + 1 ;

    type crs ;

    crs.m_memory.allocate( offset_count , label );
    crs.m_column.m_map = crs.m_memory.ptr_on_device() + row_count + 1 ;
    crs.m_row_count    = row_count ;
    crs.m_entry_count  = total_count ;

    // If same memory space then a view:
    mirror_type tmp = create_mirror::create( crs.m_memory , offset_count );

    tmp[0] = 0 ;
    total_count = 0 ;
    offset_count = row_count + 1 ;
    for ( size_t i = 0 ; i < row_count ; ++i ) {
      const std::vector< SizeType > & column = input[i] ;
      for ( size_t j = 0 ; j < column.size() ; ++j , ++offset_count ) {
        tmp[offset_count] = column[j] ;
      }
      tmp[i+1] = total_count += column.size();
    }

    DeepCopy< memory_type , mirror_type >::run( crs.m_memory , tmp , offset_count );

    return crs ;
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_CRSMAP_CREATE_HPP */


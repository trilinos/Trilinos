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

#ifndef KOKKOS_HOST_MEMORYMANAGER_HPP
#define KOKKOS_HOST_MEMORYMANAGER_HPP

#include <iosfwd>
#include <typeinfo>
#include <string>

#include <impl/KokkosArray_forward.hpp>
#include <impl/KokkosArray_MemoryView.hpp>
#include <impl/KokkosArray_ViewTracker.hpp>

/*--------------------------------------------------------------------------*/

namespace KokkosArray {
namespace Impl {

/** \brief  Memory management on the host for devices */

template<>
class MemoryManager< Host > {
public:

  static int m_memory_view_tracking ;

public:

  typedef ViewTracker view_tracker ;

  static void deallocate( void * );

  static void * allocate( const std::string & label ,
                          const std::type_info & type ,
                          const size_t member_size ,
                          const size_t member_count );

  static
  void init( ViewTracker & tracker )
    { tracker.next = 0 ; }

  static
  void clear( ViewTracker & tracker , void * ptr_on_device )
    {
      if ( tracker.remove_and_query_is_last() ) {
        deallocate( ptr_on_device );
      }
    }

  static
  void assign( ViewTracker & lhs , const ViewTracker & rhs )
  {
    if ( m_memory_view_tracking ) {
      lhs.insert( rhs );
    }
  }

  static void print_memory_view( std::ostream & );

  /*--------------------------------*/

  static void disable_memory_view_tracking();
  static void enable_memory_view_tracking();

  /*--------------------------------*/

  template < typename ValueType >
  static 
  size_t preferred_alignment( size_t parallel_length )
  {
    const size_t page = Host::detect_memory_page_size();
    if ( 0 == page % sizeof(ValueType) ) {
      const size_t align = page / sizeof(ValueType);
      const size_t rem   = parallel_length % align ;
      if ( rem ) parallel_length += align - rem ;
    }
    return parallel_length ;
  }
};

//----------------------------------------------------------------------------

template< typename ValueType >
struct Factory< MemoryView< ValueType, Host > ,
                MemoryView< ValueType, Host > >
{
  static inline
  void deep_copy( const MemoryView< ValueType , Host > & dst ,
                  const MemoryView< ValueType , Host > & src ,
                  const size_t count )
  {
    ValueType * d = dst.ptr_on_device();
    const ValueType * s = src.ptr_on_device();
    if ( d != s ) {
      const ValueType * const s_end = s + count ;
      while ( s_end != s ) { *d++ = *s++ ; }
    }
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#include <KokkosArray_Host_macros.hpp>
#include <impl/KokkosArray_MemoryView_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>

#endif /* #define KOKKOS_HOST_MEMORYMANAGER_HPP */


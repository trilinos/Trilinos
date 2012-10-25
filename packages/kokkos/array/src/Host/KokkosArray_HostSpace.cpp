/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <memory.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

#include <KokkosArray_HostSpace.hpp>
#include <Host/KokkosArray_Host_Internal.hpp>
#include <impl/KokkosArray_Error.hpp>
#include <impl/KokkosArray_MemoryTracking.hpp>

/*--------------------------------------------------------------------------*/

namespace KokkosArray {
namespace {

Impl::MemoryTracking & host_space_singleton()
{
  static Impl::MemoryTracking self("KokkosArray::HostSpace");
  return self ;
}

} // namespace <blank>
} // namespade KokkosArray

/*--------------------------------------------------------------------------*/

namespace KokkosArray {

void * HostSpace::allocate(
  const std::string    & label ,
  const std::type_info & scalar_type ,
  const size_t           scalar_size ,
  const size_t           scalar_count )
{
  assert_master_thread( "KokkosArray::HostSpace::allocate" );

  void * ptr = 0 ;

  if ( 0 < scalar_size * scalar_count ) {

    ptr = malloc( scalar_size * scalar_count );

    if ( 0 == ptr ) {
      std::ostringstream msg ;
      msg << "KokkosArray::Impl::HostSpace::allocate( "
          << label
          << " , " << scalar_type.name()
          << " , " << scalar_size
          << " , " << scalar_count
          << " ) FAILED memory allocation" ;
      KokkosArray::Impl::throw_runtime_exception( msg.str() );
    }

    host_space_singleton()
      .track( ptr, & scalar_type, scalar_size, scalar_count, label );
  }

  return ptr ;
}

void HostSpace::increment( const void * ptr )
{
  if ( 0 != ptr && Impl::HostInternal::singleton().is_master_thread() ) {
    host_space_singleton().increment( ptr );
  }
}

void HostSpace::decrement( const void * ptr )
{
  if ( 0 != ptr && Impl::HostInternal::singleton().is_master_thread() ) {

    void * ptr_alloc = host_space_singleton().decrement( ptr );

    if ( 0 != ptr_alloc ) {
      free( ptr_alloc );
    }
  }
}

void HostSpace::print_memory_view( std::ostream & o )
{
  host_space_singleton().print( o , std::string("  ") );
}

std::string HostSpace::query_label( const void * p )
{
  const Impl::MemoryTracking::Info info = 
    host_space_singleton().query( p );

  return info.label ;
}

size_t HostSpace::preferred_alignment(
  size_t scalar_size , size_t scalar_count )
{
  const size_t alignment = Host::detect_cache_line_size();

  // If the array is larger than the cache line
  // then align the count on cache line boundary.

  if ( alignment < scalar_size * scalar_count &&
       0 == alignment % scalar_size ) {
    const size_t align = alignment / scalar_size ;
    const size_t rem   = scalar_count % align ;
    if ( rem ) scalar_count += align - rem ;
  }
  return scalar_count ;
}


DeepCopy<HostSpace,HostSpace>
  ::DeepCopy( void * dst , const void * src , size_t n )
{
  memcpy( dst , src , n );
}


} // namespace KokkosArray

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace KokkosArray {

size_t HostSpace::detect_cache_line_size()
{ return Impl::HostInternal::singleton().m_cache_line_size ; }

size_t HostSpace::detect_memory_page_size()
{ return Impl::HostInternal::singleton().m_page_size ; }

void HostSpace::assert_master_thread( const char * const name )
{
  if ( ! Impl::HostInternal::singleton().is_master_thread() ) {
    std::string msg ;
    msg.append( "KokkosArray::HostSpace::assert_master_thread( " );
    msg.append( name );
    msg.append( " ) FAILED " );
    KokkosArray::Impl::throw_runtime_exception( msg );
  }
}

} // namespace KokkosArray

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/


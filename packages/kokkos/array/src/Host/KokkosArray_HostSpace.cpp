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
#include <stddef.h>
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
namespace Impl {

void * host_allocate_not_thread_safe(
  const std::string    & label ,
  const std::type_info & scalar_type ,
  const size_t           scalar_size ,
  const size_t           scalar_count )
{
  void * ptr = 0 ;

  if ( 0 < scalar_size && 0 < scalar_count ) {
    void * ptr_alloc = 0 ;
    size_t count_alloc = scalar_count ;

#if defined( __INTEL_COMPILER )

    ptr = ptr_alloc = _mm_malloc( scalar_size * count_alloc , HostSpace::MEMORY_ALIGNMENT );
   
#elif ( defined( _POSIX_C_SOURCE ) && _POSIX_C_SOURCE >= 200112L ) || \
      ( defined( _XOPEN_SOURCE )   && _XOPEN_SOURCE   >= 600 )

    posix_memalign( & ptr_alloc , HostSpace::MEMORY_ALIGNMENT , scalar_size * count_alloc );
    ptr = ptr_alloc ;

#else

    // Over-allocate to guarantee enough aligned space.

    count_alloc += ( HostSpace::MEMORY_ALIGNMENT + scalar_size - 1 ) / scalar_size ;

    ptr_alloc = malloc( scalar_size * count_alloc );

    ptr = static_cast<unsigned char *>(ptr_alloc) + 
          ( HostSpace::MEMORY_ALIGNMENT - reinterpret_cast<ptrdiff_t>(ptr_alloc) % HostSpace::MEMORY_ALIGNMENT );

#endif

    if ( ptr_alloc && ptr_alloc <= ptr &&
         0 == ( reinterpret_cast<ptrdiff_t>(ptr) % HostSpace::MEMORY_ALIGNMENT ) ) {
      host_space_singleton()
        .track( ptr_alloc, & scalar_type, scalar_size, count_alloc, label );
    }
    else {
      std::ostringstream msg ;
      msg << "KokkosArray::Impl::host_allocate_not_thread_safe( "
          << label
          << " , " << scalar_type.name()
          << " , " << scalar_size
          << " , " << scalar_count
          << " ) FAILED aligned memory allocation" ;
      KokkosArray::Impl::throw_runtime_exception( msg.str() );
    }
  }

  return ptr ;
}

void host_decrement_not_thread_safe( const void * ptr )
{
  if ( 0 != ptr ) {

    void * ptr_alloc = host_space_singleton().decrement( ptr );

    if ( 0 != ptr_alloc ) {
#if defined( __INTEL_COMPILER )
       _mm_free( ptr_alloc );
#else
       free( ptr_alloc );
#endif
    }
  }
}

}
}

namespace KokkosArray {

void * HostSpace::allocate(
  const std::string    & label ,
  const std::type_info & scalar_type ,
  const size_t           scalar_size ,
  const size_t           scalar_count )
{
  assert_master_thread( "KokkosArray::HostSpace::allocate" );

  void * const ptr =
    Impl::host_allocate_not_thread_safe( label , scalar_type , scalar_size , scalar_count );

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
  if ( Impl::HostInternal::singleton().is_master_thread() ) {
    Impl::host_decrement_not_thread_safe( ptr );
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
  // Padding 'large' array dimensions to be memory aligned
  // where 'large' greater than 4 * cacheline-size

  const size_t align = 0 == MEMORY_ALIGNMENT % scalar_size
                     ? MEMORY_ALIGNMENT / scalar_size : 0 ;

  const size_t threshold = align * 4 ;

  if ( align && threshold < scalar_count && scalar_count % align ) {
    scalar_count += align - scalar_count % align ;
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


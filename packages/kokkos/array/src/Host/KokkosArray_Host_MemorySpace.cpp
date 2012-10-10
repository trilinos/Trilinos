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

#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <sstream>

#include <KokkosArray_Host.hpp>
#include <Host/KokkosArray_Host_MemorySpace.hpp>
#include <Host/KokkosArray_Host_Internal.hpp>
#include <impl/KokkosArray_MemoryTracking.hpp>

/*--------------------------------------------------------------------------*/

namespace KokkosArray {
namespace Impl {

namespace {

class HostMemoryImpl {
public:
  MemoryTracking m_allocations ;

  HostMemoryImpl();
  ~HostMemoryImpl();

  static HostMemoryImpl & singleton();
};

HostMemoryImpl::HostMemoryImpl()
  : m_allocations()
{}

HostMemoryImpl & HostMemoryImpl::singleton()
{
  static HostMemoryImpl self ;
  return self ;
}

HostMemoryImpl::~HostMemoryImpl()
{
  if ( ! m_allocations.empty() ) {
    std::cerr << "KokkosArray::Host memory leaks:" << std::endl ;
    m_allocations.print( std::cerr , std::string("  ") );
  }
}

}

/*--------------------------------------------------------------------------*/

void * HostMemorySpace::allocate(
  const std::string    & label ,
  const std::type_info & scalar_type ,
  const size_t           scalar_size ,
  const size_t           scalar_count )
{
  HostMemoryImpl & s = HostMemoryImpl::singleton();

  if ( ! HostInternal::singleton().is_master_thread() ) {
    std::ostringstream msg ;
    msg << "KokkosArray::Impl::HostMemorySpace::allocate FAILED "
        << "can only be called on the master thread." ;
    throw std::runtime_error( msg.str() );
  }

  void * ptr = 0 ;

  if ( 0 < scalar_size * scalar_count ) {

    ptr = malloc( scalar_size * scalar_count );

    if ( 0 == ptr ) {
      std::ostringstream msg ;
      msg << "KokkosArray::Impl::HostMemorySpace::allocate( "
          << label
          << " , " << scalar_type.name()
          << " , " << scalar_size
          << " , " << scalar_count
          << " ) FAILED memory allocation" ;
      throw std::runtime_error( msg.str() );
    }

    s.m_allocations.track( ptr, & scalar_type, scalar_size, scalar_count, label );
  }

  return ptr ;
}

void HostMemorySpace::increment( const void * ptr )
{
  if ( 0 != ptr && HostInternal::singleton().is_master_thread() ) {
    HostMemoryImpl & s = HostMemoryImpl::singleton();

    s.m_allocations.increment( ptr );
  }
}

void HostMemorySpace::decrement( const void * ptr )
{
  if ( 0 != ptr && HostInternal::singleton().is_master_thread() ) {
    HostMemoryImpl & s = HostMemoryImpl::singleton();

    void * ptr_alloc = s.m_allocations.decrement( ptr );

    if ( 0 != ptr_alloc ) {
      free( ptr_alloc );
    }
  }
}

void HostMemorySpace::print_memory_view( std::ostream & o )
{
  HostMemoryImpl & s = HostMemoryImpl::singleton();

  s.m_allocations.print( o , std::string("  ") );
}


size_t HostMemorySpace::preferred_alignment(
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

/*--------------------------------------------------------------------------*/

} // namespace Impl
} // namespace KokkosArray


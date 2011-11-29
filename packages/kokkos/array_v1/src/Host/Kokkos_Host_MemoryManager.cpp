/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <sstream>

#include <Kokkos_Host.hpp>
#include <impl/Kokkos_MemoryInfo.hpp>
#include <Host/Kokkos_Host_Internal.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

namespace {

class HostMemoryImpl {
public:
  Impl::MemoryInfoSet m_allocations ;
  const size_t        m_page_size ;

  HostMemoryImpl();
  ~HostMemoryImpl();

  static HostMemoryImpl & singleton();
};

HostMemoryImpl::HostMemoryImpl()
  : m_allocations()
  , m_page_size( host_internal_page_size() )
{}

HostMemoryImpl & HostMemoryImpl::singleton()
{
  static HostMemoryImpl self ;
  return self ;
}

HostMemoryImpl::~HostMemoryImpl()
{
  if ( ! m_allocations.empty() ) {
    std::cerr << "Kokkos::Host memory leaks:" << std::endl ;
    m_allocations.print( std::cerr );
  }
}

}

/*--------------------------------------------------------------------------*/

void * MemoryManager< Host >::allocate(
  const std::string    & label ,
  const std::type_info & type ,
  const size_t member_size ,
  const size_t member_count )
{
  HostMemoryImpl & s = HostMemoryImpl::singleton();

  Impl::MemoryInfo tmp ;

  tmp.m_type  = & type ;
  tmp.m_label = label ;
  tmp.m_size  = member_size ;
  tmp.m_count = member_count ;
  tmp.m_ptr   = malloc( member_size * member_count );

  const bool ok_alloc  = 0 != tmp.m_ptr ;
  const bool ok_insert = ok_alloc && s.m_allocations.insert( tmp );

  if ( ! ok_alloc || ! ok_insert ) {
    std::ostringstream msg ;
    msg << "Kokkos::Impl::MemoryManager<Host>::allocate( " << label
        << " , " << type.name()
        << " , " << member_size
        << " , " << member_count
        << " ) FAILED " ;
    if ( ok_alloc ) { msg << "memory allocation" ; }
    else            { msg << "with internal error" ; }
    throw std::runtime_error( msg.str() );
  }

  return tmp.m_ptr ;
}

void MemoryManager< Host >::deallocate( void * ptr )
{
  HostMemoryImpl & s = HostMemoryImpl::singleton();

  if ( ! s.m_allocations.erase( ptr ) ) {
    std::ostringstream msg ;
    msg << "Kokkos::Impl::MemoryManager<Host>::deallocate( " << ptr
        << " ) FAILED memory allocated by this device" ;
    throw std::runtime_error( msg.str() );
  }

  free( ptr );
}

void MemoryManager< Host >::print_memory_view( std::ostream & o )
{
  HostMemoryImpl & s = HostMemoryImpl::singleton();

  s.m_allocations.print( o );
}

/*--------------------------------------------------------------------------*/

int MemoryManager< Host >::m_memory_view_tracking = true ;

void MemoryManager< Host >::disable_memory_view_tracking()
{
  if ( ! m_memory_view_tracking ) {
    std::string msg ;
    msg.append( "Kokkos::Impl::HostMemory::disable_memory_view_tracking ");
    msg.append( " FAILED memory_view_tracking already disabled" );
    throw std::runtime_error( msg );
  }
  m_memory_view_tracking = false ;
}

void MemoryManager< Host >::enable_memory_view_tracking()
{
  if ( m_memory_view_tracking ) {
    std::string msg ;
    msg.append( "Kokkos::Impl::HostMemory::enable_memory_view_tracking ");
    msg.append( " FAILED memory_view_tracking already enabled" );
    throw std::runtime_error( msg );
  }
  m_memory_view_tracking = true ;
}

/*--------------------------------------------------------------------------*/

size_t MemoryManager< Host >::detect_memory_page_size()
{ return HostMemoryImpl::singleton().m_page_size ; }

/*--------------------------------------------------------------------------*/

} // namespace Impl
} // namespace Kokkos


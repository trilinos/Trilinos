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
#include <map>
#include <ostream>
#include <Kokkos_HostTPI.hpp>

namespace Kokkos {
namespace {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class HostTPIImpl {
public:
  std::map<void*,std::string> m_allocations ;

  // Appropriate cached device information

  HostTPIImpl();

  static HostTPIImpl & singleton();
};

HostTPIImpl::HostTPIImpl()
  : m_allocations()
{
  // Appropriate device queries

}

HostTPIImpl & HostTPIImpl::singleton()
{
  static HostTPIImpl * impl = NULL ;
  if ( impl == NULL ) { impl = new HostTPIImpl(); }
  return *impl ;
}

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void * HostTPI::allocate_memory( size_type member_size ,
                                    size_type member_count ,
                                    const std::string & label )
{
  void * ptr_on_device = NULL ;

  ptr_on_device = calloc( member_size , member_count );

  HostTPIImpl::singleton().m_allocations[ ptr_on_device ] = label ;

  return ptr_on_device ;
}

void HostTPI::deallocate_memory( void * ptr_on_device )
{
  free( ptr_on_device );

  HostTPIImpl::singleton().m_allocations.erase( ptr_on_device );
}

void HostTPI::print_allocations( std::ostream & s ) const
{
  HostTPIImpl & impl = HostTPIImpl::singleton();

  std::map<void*,std::string>::const_iterator i = impl.m_allocations.begin();
  std::map<void*,std::string>::const_iterator end = impl.m_allocations.end();

  for ( ; i != end ; ++i ) {
    s << i->second << std::endl ;
  }
}


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Kokkos




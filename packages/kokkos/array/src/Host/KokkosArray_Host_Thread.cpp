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

#include <Host/KokkosArray_Host_Thread.hpp>

#include <limits>
#include <stdexcept>
#include <sstream>
#include <iostream>

namespace KokkosArray {
namespace Impl {

class HostThreadSentinel {
public:
  HostThreadSentinel();
  ~HostThreadSentinel();
  static void singleton();
};

void HostThread::set_thread( const unsigned global_rank , HostThread * t )
{
  HostThreadSentinel::singleton();

  const bool ok_rank = global_rank < max_thread_count ;
  const bool ok_zero = ok_rank && ( 0 == m_thread[ global_rank ] );

  if ( ok_rank && ok_zero ) {
    m_thread[ global_rank ] = t ;
  }
  else {
    std::ostringstream msg ;
    msg << "KokkosArray::Impl::HostThread::set_thread( "
        << global_rank << " , ... ) ERROR: " ;
    if ( ! ok_rank ) { msg << " OUT OF BOUNDS" ; }
    else if ( ! ok_zero ) { msg << " ALREADY SET" ; }
    throw std::runtime_error( msg.str() );
  }
}

void HostThread::clear_thread( const unsigned global_rank )
{
  HostThreadSentinel::singleton();

  const bool ok_rank = global_rank < max_thread_count ;

  if ( ok_rank ) {
    m_thread[ global_rank ] = 0 ;
  }
  else {
    std::ostringstream msg ;
    msg << "KokkosArray::Impl::HostThread::clear_thread( "
        << global_rank << " , ... ) ERROR:  OUT OF BOUNDS" ;
    throw std::runtime_error( msg.str() );
  }
}

HostThread::HostThread()
{
  m_fan_count    = 0 ;
  m_thread_rank  = std::numeric_limits<unsigned>::max();
  m_thread_count = 0 ;
  m_gang_rank    = std::numeric_limits<unsigned>::max();
  m_gang_count   = 0 ;
  m_worker_rank  = std::numeric_limits<unsigned>::max();
  m_worker_count = 0 ;
  m_reduce       = 0 ;
  m_state        = ThreadActive ;

  for ( unsigned i = 0 ; i < max_fan_count ; ++i ) { m_fan[i] = 0 ; }
}

HostThread::~HostThread()
{
  m_fan_count    = 0 ;
  m_thread_rank  = std::numeric_limits<unsigned>::max();
  m_thread_count = 0 ;
  m_gang_rank    = std::numeric_limits<unsigned>::max();
  m_gang_count   = 0 ;
  m_worker_rank  = std::numeric_limits<unsigned>::max();
  m_worker_count = 0 ;
  m_reduce       = 0 ;

  for ( unsigned i = 0 ; i < max_fan_count ; ++i ) { m_fan[i] = 0 ; }
}

HostThread * HostThread::m_thread[ HostThread::max_thread_count ];

//----------------------------------------------------------------------------

void HostThreadSentinel::singleton()
{
  static HostThreadSentinel self ;
}

HostThreadSentinel::HostThreadSentinel()
{
  for ( unsigned i = 0 ; i < HostThread::max_thread_count ; ++i ) {
    HostThread::m_thread[i] = 0 ;
  }
}

HostThreadSentinel::~HostThreadSentinel()
{
  unsigned nonzero_count = 0 ;
  for ( unsigned i = 0 ; i < HostThread::max_thread_count ; ++i ) {
    if ( 0 != HostThread::m_thread[i] ) ++nonzero_count ;
  }
  if ( nonzero_count ) {
    std::cerr << "KokkosArray::Impl::HostThread WARNING Terminating with "
              << nonzero_count
              << " non-null threads."
              << std::endl ;
  }
}

} // namespace Impl
} // namespace KokkosArray



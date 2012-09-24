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

/*--------------------------------------------------------------------------*/
/* KokkosArray interfaces */

#include <KokkosArray_Host.hpp>
#include <Host/KokkosArray_Host_Internal.hpp>

/*--------------------------------------------------------------------------*/
/* Standard 'C' libraries */
#include <stdlib.h>

/* Standard 'C++' libraries */
#include <limits>
#include <iostream>
#include <stdexcept>
#include <sstream>

/*--------------------------------------------------------------------------*/

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

HostInternal::~HostInternal()
{}

HostInternal::HostInternal()
  : m_worker_block()
  , m_node_rank( -1 )
  , m_node_count( 1 )
  , m_node_pu_count( 0 )
  , m_page_size( 0 )
  , m_cache_line_size( 64 /* default alignment */ )
  , m_thread_count( 1 )
  , m_gang_count( 1 )
  , m_worker_count( 1 )
{
  m_worker = NULL ;

  if ( ! is_master_thread() ) {
    throw std::runtime_error( std::string("KokkosArray::Impl::HostInternal FAILED : not initialized on the master thread"));
  }

  // Master thread:
  m_thread[0] = & m_master_thread ;

  m_master_thread.m_fan_count    = 0 ;
  m_master_thread.m_thread_rank  = 0 ;
  m_master_thread.m_thread_count = 1 ;
  m_master_thread.m_gang_rank    = 0 ;
  m_master_thread.m_gang_count   = 1 ;
  m_master_thread.m_worker_rank  = 0 ;
  m_master_thread.m_worker_count = 1 ;

  for ( unsigned i = 0 ; i < HostThread::max_fan_count ; ++i ) {
    m_master_thread.m_fan[i] = 0 ;
  }
}

//----------------------------------------------------------------------------

bool HostInternal::initialize_thread(
  const unsigned thread_rank ,
  HostThread & thread )
{
  const unsigned gang_rank   = thread_rank / m_worker_count ;
  const unsigned worker_rank = thread_rank % m_worker_count ;

  thread.m_thread_rank  = thread_rank ;
  thread.m_thread_count = m_thread_count ;
  thread.m_worker_rank  = worker_rank ;
  thread.m_worker_count = m_worker_count ;
  thread.m_gang_rank    = gang_rank ;
  thread.m_gang_count   = m_gang_count ;

  {
    unsigned fan_count = 0 ;

    // Intranode reduction:
    for ( unsigned n = 1 ; worker_rank + n < m_worker_count ; n <<= 1 ) {

      if ( n & worker_rank ) break ;

      HostThread * const th = m_thread[ thread_rank + n ];

      if ( 0 == th ) return false ;

      thread.m_fan[ fan_count++ ] = th ;
    }

    if ( worker_rank == 0 ) {

      // Internode reduction:
      for ( unsigned n = 1 ; gang_rank + n < m_gang_count ; n <<= 1 ) {

        if ( n & gang_rank ) break ;

        HostThread * const th = m_thread[ thread_rank + n * m_worker_count ];

        if ( 0 == th ) return false ;

        thread.m_fan[ fan_count++ ] = th ;
      }
    }

    thread.m_fan_count = fan_count ;
  }

  return true ;
}

bool HostInternal::bind_thread( const unsigned thread_rank ) const
{
  (void) thread_rank; // Prevent compiler warning for unused argument
  return true ;
}

//----------------------------------------------------------------------------

void HostInternal::finalize()
{
  verify_inactive("finalize()");

  // Release and clear worker threads:
  while ( 1 < m_thread_count ) {
    --m_thread_count ;

    if ( m_thread[ m_thread_count ] ) {
      m_master_thread.set( HostThread::ThreadInactive );

      m_thread[ m_thread_count ]->set( HostThread::ThreadTerminating );

      m_master_thread.wait( HostThread::ThreadInactive );

      // Is in the 'ThreadTerminating" state
    }
  }

  // Reset master thread:
  m_master_thread.m_fan_count    = 0 ;
  m_master_thread.m_thread_rank  = 0 ;
  m_master_thread.m_thread_count = 1 ;
  m_master_thread.m_gang_rank    = 0 ;
  m_master_thread.m_gang_count   = 1 ;
  m_master_thread.m_worker_rank  = 0 ;
  m_master_thread.m_worker_count = 1 ;
  m_master_thread.set( HostThread::ThreadActive );

  for ( unsigned i = 0 ; i < HostThread::max_fan_count ; ++i ) {
    m_master_thread.m_fan[i] = 0 ;
  }
}

//----------------------------------------------------------------------------
// Driver for each created thread

void HostInternal::driver( const size_t thread_rank )
{
  // Bind this thread to a unique processing unit
  // with all members of a gang in the same NUMA region.

  if ( bind_thread( thread_rank ) ) {

    HostThread this_thread ;

    m_thread[ thread_rank ] = & this_thread ;

    // Initialize thread ranks and fan-in relationships:

    if ( initialize_thread( thread_rank , this_thread ) ) {

      // Inform master thread that binding and initialization succeeded.
      m_master_thread.set( HostThread::ThreadActive );

      try {
        // Work loop:

        while ( HostThread::ThreadActive == this_thread.m_state ) {

          // When the work is complete the state will be Inactive or Terminate
          m_worker->execute_on_thread( this_thread );

          // If this_thread is in the Inactive state then wait for activation.
          this_thread.wait( HostThread::ThreadInactive );
        }
      }
      catch( const std::exception & x ) {
        // mfh 29 May 2012: Doesn't calling std::terminate() seem a
        // little violent?  On the other hand, C++ doesn't define how
        // to transport exceptions between threads (until C++11).
        // Since this is a worker thread, it would be hard to tell the
        // master thread what happened.
        std::cerr << "Thread " << thread_rank << " uncaught exception : "
                  << x.what() << std::endl ;
        std::terminate();
      }
      catch( ... ) {
        // mfh 29 May 2012: See note above on std::terminate().
        std::cerr << "Thread " << thread_rank << " uncaught exception"
                  << std::endl ;
        std::terminate();
      }
    }
  }

  // Notify master thread that this thread has terminated.

  m_thread[ thread_rank ] = 0 ;

  m_master_thread.set( HostThread::ThreadTerminating );
}

//----------------------------------------------------------------------------

void HostInternal::verify_inactive( const char * const method ) const
{
  if ( NULL != m_worker ) {

    std::ostringstream msg ;
    msg << "KokkosArray::Host::" << method << " FAILED: " ;

    if ( & m_worker_block == m_worker ) {
      msg << "Device is blocked" ;
    }
    else {
      msg << "Functor is executing." ;
    }
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------------

bool HostInternal::spawn_threads( const unsigned gang_count ,
                                  const unsigned worker_count )
{
  // If the process is bound to a particular node
  // then only use cores belonging to that node.
  // Otherwise use all nodes and all their cores.

  m_gang_count   = gang_count ;
  m_worker_count = worker_count ;
  m_thread_count = gang_count * worker_count ;
  m_worker       = & m_worker_block ;

  // Bind the process thread as thread_rank == 0
  bool ok_spawn_threads = bind_thread( 0 );

  // Spawn threads from last-to-first so that the
  // fan-in barrier thread relationships can be established.

  for ( unsigned rank = m_thread_count ; ok_spawn_threads && 0 < --rank ; ) {

    m_master_thread.set( HostThread::ThreadInactive );

    // Spawn thread executing the 'driver' function.
    ok_spawn_threads = spawn( rank );

    if ( ok_spawn_threads ) {

      // Thread spawned, wait for thread to activate:
      m_master_thread.wait( HostThread::ThreadInactive );

      // Check if the thread initialized and bound correctly:
      ok_spawn_threads = HostThread::ThreadActive == m_master_thread.m_state ;

      if ( ok_spawn_threads ) { // Wait for spawned thread to deactivate
        HostThread * volatile * const threads = m_thread ;
        threads[ rank ]->wait( HostThread::ThreadActive );
        // m_thread[ rank ]->wait( HostThread::ThreadActive );
      }
    }
  }

  m_worker = NULL ;

  // All threads spawned, initialize the master-thread fan-in
  ok_spawn_threads =
    ok_spawn_threads && initialize_thread( 0 , m_master_thread );

  if ( ! ok_spawn_threads ) {
    finalize();
  }

  return ok_spawn_threads ;
}

void HostInternal::initialize( const unsigned gang_count ,
                               const unsigned worker_count )
{
  const bool ok_inactive     = 1 == m_thread_count ;
  const bool ok_gang_count   = gang_count <= m_node_count ;
  const bool ok_worker_count = ( 0 == m_node_pu_count ||
                                 worker_count <= m_node_pu_count );

  // Only try to spawn threads if input is valid.

  const bool ok_spawn_threads =
    ( ok_inactive &&
      ok_gang_count && ok_worker_count &&
      1 < gang_count * worker_count )
    ? spawn_threads( gang_count , worker_count )
    : true ;

  if ( ! ok_inactive ||
       ! ok_gang_count ||
       ! ok_worker_count ||
       ! ok_spawn_threads )
  {
    std::ostringstream msg ;

    msg << "KokkosArray::Host::initialize() FAILED" ;

    if ( ! ok_inactive ) {
      msg << " : Device is already active" ;
    }
    if ( ! ok_gang_count ) {
      msg << " : gang_count(" << gang_count
          << ") exceeds detect_node_count(" << m_node_count
          << ")" ;
    }
    if ( ! ok_worker_count ) {
      msg << " : worker_count(" << worker_count
          << ") exceeds detect_node_pu_count(" << m_node_pu_count
          << ")" ;
    }
    if ( ! ok_spawn_threads ) {
      msg << " : Spawning or cpu-binding the threads" ;
    }

    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------------

void HostInternal::activate()
{
  for ( unsigned i = m_thread_count ; 1 < i ; ) {
    m_thread[--i]->set( HostThread::ThreadActive );
  }
}

inline
void HostInternal::execute( const HostThreadWorker & worker )
{
  verify_inactive("execute(...)");

  m_worker = & worker ;

  activate();

  // Will finalize with a barrier.
  worker.execute_on_thread( m_master_thread );

  m_worker = NULL ;
}

void HostThreadWorker::execute( const HostThreadWorker & worker )
{ HostInternal::singleton().execute( worker ); }

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

void Host::finalize()
{ Impl::HostInternal::singleton().finalize(); }

void Host::initialize( const Host::size_type gang_count ,
                       const Host::size_type worker_count )
{
  Impl::HostInternal::singleton().initialize( gang_count , worker_count );
}

Host::size_type Host::detect_node_count()
{ return Impl::HostInternal::singleton().m_node_count ; }

Host::size_type Host::detect_node_core_count()
{ return Impl::HostInternal::singleton().m_node_pu_count ; }

Host::size_type Host::detect_cache_line_size()
{ return Impl::HostInternal::singleton().m_cache_line_size ; }

Host::size_type Host::detect_memory_page_size()
{ return Impl::HostInternal::singleton().m_page_size ; }

/*--------------------------------------------------------------------------*/

} // namespace KokkosArray


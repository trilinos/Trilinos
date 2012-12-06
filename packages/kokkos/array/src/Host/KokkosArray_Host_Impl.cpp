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
#include <impl/KokkosArray_Error.hpp>

/*--------------------------------------------------------------------------*/
/* Standard 'C' libraries */
#include <stdlib.h>

/* Standard 'C++' libraries */
#include <limits>
#include <iostream>
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
  m_work_chunk   = 0 ;
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
  m_work_chunk   = 0 ;
  m_reduce       = 0 ;
  m_state        = ThreadActive ;

  for ( unsigned i = 0 ; i < max_fan_count ; ++i ) { m_fan[i] = 0 ; }
}

std::pair< Host::size_type , Host::size_type >
HostThread::work_range( const size_type work_count ) const
{
  const size_type chunk_count =
    ( work_count + m_work_chunk - 1 ) / m_work_chunk ;

  const size_type work_per_thread =
    m_work_chunk * (( chunk_count + m_thread_count - 1 ) / m_thread_count );

  const size_type work_begin =
    std::min( m_thread_rank * work_per_thread , work_count );

  const size_type work_end =
    std::min( work_begin + work_per_thread , work_count );

  return std::pair<size_type,size_type>( work_begin , work_end );
}

void HostThread::barrier()
{
  // The 'wait' function repeatedly polls the 'thread' state
  // which may reside in a different NUMA region.
  // Thus the fan is intra-node followed by inter-node
  // to minimize inter-node memory access.

  for ( unsigned i = 0 ; i < m_fan_count ; ++i ) {
    // Wait until thread enters the 'Rendezvous' state
    m_fan[i]->wait( HostThread::ThreadActive );
  }

  if ( m_thread_rank ) {
    set(  HostThread::ThreadRendezvous );
    wait( HostThread::ThreadRendezvous );
  }

  for ( unsigned i = m_fan_count ; 0 < i ; ) {
    m_fan[--i]->set( HostThread::ThreadActive );
  }
}

/** \brief  This thread waits for each fan-in thread to deactivate.  */
void HostThreadWorker::fanin_deactivation( HostThread & thread ) const
{
  // The 'wait' function repeatedly polls the 'thread' state
  // which may reside in a different NUMA region.
  // Thus the fan is intra-node followed by inter-node
  // to minimize inter-node memory access.

  for ( unsigned i = 0 ; i < thread.m_fan_count ; ++i ) {
    thread.m_fan[i]->wait( HostThread::ThreadActive );
  }
}

void HostInternal::activate_threads()
{
  // Activate threads to call 'm_worker.execute_on_thread'
  for ( unsigned i = m_thread_count ; 1 < i ; ) {
    m_thread[--i]->set( HostThread::ThreadActive );
  }
}

inline
void HostInternal::execute( const HostThreadWorker & worker )
{
  verify_inactive("execute(...)");

  // Worker threads are verified to be in the ThreadInactive state.

  m_worker = & worker ;

  activate_threads();

  // Execute on the master thread,
  worker.execute_on_thread( m_master_thread );

  // Wait for fanin/fanout threads to self-deactivate.
  worker.fanin_deactivation( m_master_thread );

  // Worker threads are returned to the ThreadInactive state.
  m_worker = NULL ;
}

void HostThreadWorker::execute() const
{ HostInternal::singleton().execute( *this ); }

void HostInternal::execute_serial( const HostThreadWorker & worker )
{
  verify_inactive("execute_serial(...)");

  // Worker threads are verified to be in the ThreadInactive state.

  m_worker = & worker ;

  for ( unsigned i = m_thread_count ; 1 < i ; ) {
    --i ;
    m_thread[i]->set(  HostThread::ThreadActive );
    m_thread[i]->wait( HostThread::ThreadActive );
  }

  // Execute on the master thread,
  worker.execute_on_thread( m_master_thread );

  // Worker threads are returned to the ThreadInactive state.
  m_worker = NULL ;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

HostInternal::~HostInternal() {}

HostInternal::HostInternal()
  : m_worker_block()
  , m_node_rank( -1 )
  , m_node_count( 1 )
  , m_node_pu_count( 0 )
  , m_page_size( 0 )
  , m_cache_line_size( 64 /* default guess */ )
  , m_thread_count( 1 )
  , m_gang_count( 1 )
  , m_worker_count( 1 )
  , m_work_chunk( m_cache_line_size / sizeof(void*) )
  , m_reduce_scratch_size( 0 )
  , m_reduce_scratch( 0 )
{
  m_worker = NULL ;

  if ( ! is_master_thread() ) {
    KokkosArray::Impl::throw_runtime_exception( std::string("KokkosArray::Impl::HostInternal FAILED : not initialized on the master thread") );
  }

  // Master thread:
  m_thread[0] = & m_master_thread ;

  m_master_thread.m_fan_count    = 0 ;
  m_master_thread.m_thread_rank  = 0 ;
  m_master_thread.m_thread_count = m_thread_count ;
  m_master_thread.m_gang_rank    = 0 ;
  m_master_thread.m_gang_count   = m_gang_count ;
  m_master_thread.m_worker_rank  = 0 ;
  m_master_thread.m_worker_count = m_worker_count ;
  m_master_thread.m_work_chunk   = m_work_chunk ;

  for ( unsigned i = 0 ; i < HostThread::max_fan_count ; ++i ) {
    m_master_thread.m_fan[i] = 0 ;
  }
}

//----------------------------------------------------------------------------

inline
void HostInternal::resize_reduce_thread( HostThread & thread ) const
{
  if ( thread.m_reduce ) {
    free( thread.m_reduce );
    thread.m_reduce = 0 ;
  }

  if ( m_reduce_scratch_size ) {
    thread.m_reduce = malloc( m_reduce_scratch_size );

    // Guaranteed multiple of 'unsigned'

    unsigned * ptr = (unsigned *)( thread.m_reduce );
    unsigned * const end = ptr + m_reduce_scratch_size / sizeof(unsigned );
  
    // touch on this thread
    while ( ptr < end ) *ptr++ = 0 ;
  }
}

struct HostWorkerResizeReduce : public HostThreadWorker
{
  void execute_on_thread( HostThread & ) const ;

  ~HostWorkerResizeReduce() {}
  HostWorkerResizeReduce() {}
};

void HostWorkerResizeReduce::execute_on_thread( HostThread & thread ) const
{ HostInternal::singleton().resize_reduce_thread( thread ); }


inline
void * HostInternal::reduce_scratch() const
{ return m_reduce_scratch ; }

inline
void HostInternal::resize_reduce_scratch( unsigned size )
{
  if ( 0 == size || m_reduce_scratch_size < size ) {
    // Round up to cache line size:
    const unsigned rem = size % m_cache_line_size ;

    if ( rem ) size += m_cache_line_size - rem ;
  
    m_reduce_scratch_size = size ;

    const HostWorkerResizeReduce work ;

    execute_serial( work );

    m_reduce_scratch = m_master_thread.m_reduce ;
  }
}

void host_resize_scratch_reduce( unsigned size )
{ HostInternal::singleton().resize_reduce_scratch( size ); }

void * host_scratch_reduce()
{ return HostInternal::singleton().reduce_scratch(); }

//----------------------------------------------------------------------------

bool HostInternal::initialize_thread(
  const unsigned thread_rank ,
  HostThread & thread )
{
  const unsigned gang_rank   = thread_rank / m_worker_count ;
  const unsigned worker_rank = thread_rank % m_worker_count ;

  thread.m_thread_rank  = thread_rank ;
  thread.m_thread_count = m_thread_count ;
  thread.m_gang_rank    = gang_rank ;
  thread.m_gang_count   = m_gang_count ;
  thread.m_worker_rank  = worker_rank ;
  thread.m_worker_count = m_worker_count ;
  thread.m_work_chunk   = m_work_chunk ;

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

  resize_reduce_scratch( 0 );

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

          // Perform the work:
          m_worker->execute_on_thread( this_thread );

          // Wait for fanin threads to self-deactivate:
          m_worker->fanin_deactivation( this_thread );

          // Deactivate this thread:
          this_thread.set(  HostThread::ThreadInactive );

          // Wait to be activated or terminated:
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
    KokkosArray::Impl::throw_runtime_exception( msg.str() );
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
    if ( ! ok_gang_count || ! ok_worker_count ) {
      msg << " : request for threads ( "
          << gang_count << " x " << worker_count
          << " ) exceeds detected capacity ( "
          << m_node_count << " x " << m_node_pu_count
          << " )" ;
    }
    if ( ! ok_spawn_threads ) {
      msg << " : Spawning or cpu-binding the threads" ;
    }

    KokkosArray::Impl::throw_runtime_exception( msg.str() );
  }

  resize_reduce_scratch( 1024 );
}

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


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
#include <Host/KokkosArray_hwloc.hpp>
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
namespace {

class HostWorkerBlock : public HostThreadWorker {
public:

  void execute_on_thread( HostThread & ) const 
  {
    host_thread_lock();
    host_thread_unlock();
  }

  HostWorkerBlock()  {}
  ~HostWorkerBlock() {}
};

const HostWorkerBlock worker_block ;

//----------------------------------------------------------------------------

inline
void thread_mapping( const unsigned gang_rank ,
                     const unsigned gang_count,
                     const unsigned worker_rank ,
                     const unsigned worker_count ,
                     unsigned coordinate[] )
{
  unsigned capacity[ hwloc::max_depth ];

  hwloc::get_thread_capacity( capacity );

  for ( unsigned i = 0 ; i < hwloc::max_depth ; ++i ) coordinate[i] = 0 ;

  { // Assign gang to resource:
    const unsigned bin  = gang_count / capacity[0] ;
    const unsigned bin1 = bin + 1 ;
    const unsigned k    = capacity[0] * bin1 - gang_count ;
    const unsigned part = k * bin ;

    if ( gang_rank < part ) {
      coordinate[0] = gang_rank / bin ;
    }
    else {
      coordinate[0] = k + ( gang_rank - part ) / bin1 ;
    }
  }

  { // Assign workers to resources:
    unsigned n = worker_count ;
    unsigned r = worker_rank ;

    for ( unsigned i = 1 ; i < hwloc::max_depth &&
                           0 < capacity[i]  ; ++i ) {
      // n = k * bin + ( capacity[i] - k ) * ( bin + 1 )
      const unsigned bin  = n / capacity[i] ;
      const unsigned bin1 = bin + 1 ;
      const unsigned k    = capacity[i] * bin1 - n ;
      const unsigned part = k * bin ;

      if ( r < part ) {
        coordinate[i]  = r / bin ;
        r = r - bin * coordinate[i] ;
        n = bin ;
      }
      else {
        const unsigned r1 = r - part ;
        const unsigned c1 = r1 / bin1 ;

        coordinate[i]  = c1 + k ;
        r = r1 - c1 * bin1 ;
        n = bin1 ;
      }
    }
  }
}

} // namespace

//----------------------------------------------------------------------------

void host_barrier( HostThread & thread )
{
  // The 'wait' function repeatedly polls the 'thread' state
  // which may reside in a different NUMA region.
  // Thus the fan is intra-node followed by inter-node
  // to minimize inter-node memory access.

  for ( unsigned i = 0 ; i < thread.fan_count() ; ++i ) {
    // Wait until thread enters the 'Rendezvous' state
    host_thread_wait( & thread.fan(i).m_state , HostThread::ThreadActive );
  }

  if ( thread.rank() ) {
    thread.m_state = HostThread::ThreadRendezvous ;
    host_thread_wait( & thread.m_state , HostThread::ThreadRendezvous );
  }

  for ( unsigned i = thread.fan_count() ; 0 < i ; ) {
    thread.fan(--i).m_state = HostThread::ThreadActive ;
  }
}

//----------------------------------------------------------------------------

/** \brief  This thread waits for each fan-in thread to deactivate.  */
void HostThreadWorker::fanin_deactivation( HostThread & thread ) const
{
  // The 'wait' function repeatedly polls the 'thread' state
  // which may reside in a different NUMA region.
  // Thus the fan is intra-node followed by inter-node
  // to minimize inter-node memory access.

  for ( unsigned i = 0 ; i < thread.fan_count() ; ++i ) {
    host_thread_wait( & thread.fan(i).m_state , HostThread::ThreadActive );
  }
}

//----------------------------------------------------------------------------

void HostInternal::activate_threads()
{
  // Activate threads to call 'm_worker.execute_on_thread'
  for ( unsigned i = m_thread_count ; 1 < i ; ) {
    HostThread::get_thread(--i)->m_state = HostThread::ThreadActive ;
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
    HostThread & thread = * HostThread::get_thread(i);

    thread.m_state = HostThread::ThreadActive ;
    host_thread_wait( & thread.m_state , HostThread::ThreadActive );
  }

  // Execute on the master thread,
  worker.execute_on_thread( m_master_thread );

  // Worker threads are returned to the ThreadInactive state.
  m_worker = NULL ;
}

void HostThreadWorker::execute_serial() const
{ HostInternal::singleton().execute_serial( *this ); }

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

HostInternal::~HostInternal()
{
  HostThread::clear_thread(0);
}

HostInternal::HostInternal()
  : m_thread_count( 1 )
  , m_gang_count( 1 )
  , m_worker_count( 1 )
{
  m_worker = NULL ;

  if ( ! is_master_thread() ) {
    KokkosArray::Impl::throw_runtime_exception( std::string("KokkosArray::Impl::HostInternal FAILED : not initialized on the master thread") );
  }

  // Initialize thread pool with a single master thread:
  HostThread::set_thread( 0 , & m_master_thread );
  HostThread::set_thread_relationships();
}

HostInternal & HostInternal::singleton()
{
  static HostInternal self ; return self ;
}

void HostInternal::print_configuration( std::ostream & s ) const
{
  unsigned coordinate[ hwloc::max_depth ];

  hwloc::print_thread_capacity( s );

  s << std::endl ;

  for ( unsigned r = 0 ; r < m_thread_count ; ++r ) {
    const unsigned gang_rank   = r / m_worker_count ;
    const unsigned worker_rank = r % m_worker_count ;

    thread_mapping( gang_rank , m_gang_count,
                    worker_rank , m_worker_count ,
                    coordinate );

    s << "  thread[ " << r << " : (" << gang_rank << "," << worker_rank << ") ] -> bind{"
      << " " << coordinate[0]
      << " " << coordinate[1]
      << " " << coordinate[2]
      << " }" << std::endl ;
  }
}

//----------------------------------------------------------------------------

void * HostInternal::reduce_scratch() const
{
  return m_master_thread.reduce_data();
}

void host_resize_scratch_reduce( unsigned size )
{
  HostThreadResizeReduce<Host> tmp(size);
}

void * host_scratch_reduce()
{ return HostInternal::singleton().reduce_scratch(); }

//----------------------------------------------------------------------------

bool HostInternal::bind_thread( const unsigned thread_rank ) const
{
  const unsigned gang_rank   = thread_rank / m_worker_count ;
  const unsigned worker_rank = thread_rank % m_worker_count ;

  unsigned coordinate[ hwloc::max_depth ];

  thread_mapping( gang_rank , m_gang_count,
                  worker_rank , m_worker_count ,
                  coordinate );
  return hwloc::bind_this_thread( coordinate );
}

//----------------------------------------------------------------------------

void HostInternal::finalize()
{
  verify_inactive("finalize()");

  host_resize_scratch_reduce( 0 );

  // Release and clear worker threads:
  while ( 1 < m_thread_count ) {
    --m_thread_count ;

    HostThread * const thread = HostThread::get_thread( m_thread_count );

    if ( 0 != thread ) {

      m_master_thread.m_state = HostThread::ThreadInactive ;

      thread->m_state = HostThread::ThreadTerminating ;

      host_thread_wait( & m_master_thread.m_state , HostThread::ThreadInactive );

      // Is in the 'ThreadTerminating" state
    }
  }

  // Reset master thread to the "only thread" condition:
  HostThread::set_thread_relationships();
}

//----------------------------------------------------------------------------
// Driver for each created thread

void HostInternal::driver( const size_t thread_rank )
{
  // Bind this thread to a unique processing unit
  // with all members of a gang in the same NUMA region.

//  if ( bind_thread( thread_rank ) ) {
  if ( HostInternal::bind_thread( thread_rank ) ) {

    HostThread this_thread ;

    HostThread::set_thread( thread_rank , & this_thread );

    this_thread.m_state = HostThread::ThreadActive ;

    // Inform master thread that spawning and binding succeeded.

    m_master_thread.m_state = HostThread::ThreadActive ;

    try {
      // Work loop:

      while ( HostThread::ThreadActive == this_thread.m_state ) {

        // Perform the work:
        m_worker->execute_on_thread( this_thread );

        // Wait for fanin threads to self-deactivate:
        m_worker->fanin_deactivation( this_thread );

        // Deactivate this thread:
        this_thread.m_state = HostThread::ThreadInactive ;

        // Wait to be activated or terminated:
        host_thread_wait( & this_thread.m_state , HostThread::ThreadInactive );
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

  // Notify master thread that this thread has terminated.

  HostThread::clear_thread( thread_rank );

  m_master_thread.m_state = HostThread::ThreadTerminating ;
}

//----------------------------------------------------------------------------

void HostInternal::verify_inactive( const char * const method ) const
{
  if ( NULL != m_worker ) {

    std::ostringstream msg ;
    msg << "KokkosArray::Host::" << method << " FAILED: " ;

    if ( & worker_block == m_worker ) {
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
  m_worker       = & worker_block ;

  bool ok_spawn_threads = true ;

  for ( unsigned rank = gang_count * worker_count ; ok_spawn_threads && 0 < --rank ; ) {

    m_master_thread.m_state = HostThread::ThreadInactive ;

    // Spawn thread executing the 'driver' function.
    ok_spawn_threads = spawn( rank );

    if ( ok_spawn_threads ) {

      // Thread spawned, wait for thread to activate:
      host_thread_wait( & m_master_thread.m_state , HostThread::ThreadInactive );

      // Check if the thread initialized and bound correctly:
      ok_spawn_threads = HostThread::ThreadActive == m_master_thread.m_state ;

      if ( ok_spawn_threads ) { // Wait for spawned thread to deactivate
        host_thread_wait( & HostThread::get_thread(rank)->m_state , HostThread::ThreadActive );
      }
    }
  }

  // Bind the process thread as thread_rank == 0
  if ( ok_spawn_threads ) { ok_spawn_threads = HostInternal::bind_thread( 0 ); }

  m_worker = NULL ;

  return ok_spawn_threads ;
}

void HostInternal::initialize( const unsigned gang_count ,
                               const unsigned worker_count )
{
  const bool ok_inactive = 1 == m_thread_count ;

  // Only try to spawn threads if input is valid.

  const bool ok_spawn_threads =
    ( ok_inactive && 1 < gang_count * worker_count )
    ? spawn_threads( gang_count , worker_count )
    : true ;

  if ( ok_spawn_threads ) {
    // All threads spawned, initialize the master-thread fan-in

    for ( unsigned rank = 0 ; rank < m_thread_count ; ++rank ) {
      HostThread::get_thread(rank)->m_gang_tag   = rank / m_worker_count ;
      HostThread::get_thread(rank)->m_worker_tag = rank % m_worker_count ;
    }

    HostThread::set_thread_relationships();
  }
  else {
    finalize();
  }

  if ( ! ok_inactive || ! ok_spawn_threads )
  {
    std::ostringstream msg ;

    msg << "KokkosArray::Host::initialize() FAILED" ;

    if ( ! ok_inactive ) {
      msg << " : Device is already active" ;
    }
    if ( ! ok_spawn_threads ) {
      msg << " : Spawning or cpu-binding the threads" ;
    }

    KokkosArray::Impl::throw_runtime_exception( msg.str() );
  }

  host_resize_scratch_reduce( 1024 );
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
                       const Host::size_type gang_worker_count )
{
  Impl::HostInternal::singleton().initialize( gang_count , gang_worker_count );
}

void Host::print_configuration( std::ostream & s )
{
  Impl::HostInternal::singleton().Impl::HostInternal::print_configuration(s);
}

Host::size_type Host::detect_gang_capacity()
{
  unsigned capacity[ Impl::hwloc::max_depth ];

  Impl::hwloc::get_thread_capacity( capacity );

  return 0 == capacity[0] ? 1 : capacity[0] ;
}

Host::size_type Host::detect_gang_worker_capacity()
{
  unsigned capacity[ Impl::hwloc::max_depth ];

  Impl::hwloc::get_thread_capacity( capacity );

  return 0 == capacity[0] ? 1 : (
         0 == capacity[1] ? 1 : capacity[1] * (
         0 == capacity[2] ? 1 : capacity[2] * (
         0 == capacity[3] ? 1 : capacity[3] ) ) );
}

//----------------------------------------------------------------------------

bool Host::sleep()
{
  Impl::HostInternal & h = Impl::HostInternal::singleton();

  const bool is_ready   = NULL == h.m_worker ;
        bool is_blocked = & Impl::worker_block == h.m_worker ;

  if ( is_ready ) {
    Impl::host_thread_lock();

    h.m_worker = & Impl::worker_block ;

    // Activate threads so that they will proceed from
    // spinning state to being blocked on the mutex.

    h.activate_threads();

    is_blocked = true ;
  }

  return is_blocked ;
}

bool Host::wake()
{
  Impl::HostInternal & h = Impl::HostInternal::singleton();

  const bool is_blocked = & Impl::worker_block == h.m_worker ;
        bool is_ready   = NULL == h.m_worker ;

  if ( is_blocked ) {
    Impl::host_thread_unlock();

    Impl::worker_block.fanin_deactivation( h.m_master_thread );

    h.m_worker = NULL ;

    is_ready = true ;
  }

  return is_ready ;
}

/*--------------------------------------------------------------------------*/

} // namespace KokkosArray


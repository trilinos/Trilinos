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

void * host_allocate_not_thread_safe(
  const std::string    & label ,
  const std::type_info & scalar_type ,
  const size_t           scalar_size ,
  const size_t           scalar_count );

void host_decrement_not_thread_safe( const void * ptr );

//----------------------------------------------------------------------------

void host_barrier( HostThread & thread )
{
  // The 'wait' function repeatedly polls the 'thread' state
  // which may reside in a different NUMA region.
  // Thus the fan is intra-node followed by inter-node
  // to minimize inter-node memory access.

  for ( unsigned i = 0 ; i < thread.fan_count() ; ++i ) {
    // Wait until thread enters the 'Rendezvous' state
    host_wait( & thread.fan(i).m_state , HostThread::ThreadActive );
  }

  if ( thread.rank() ) {
    thread.m_state = HostThread::ThreadRendezvous ;
    host_wait( & thread.m_state , HostThread::ThreadRendezvous );
  }

  for ( unsigned i = thread.fan_count() ; 0 < i ; ) {
    thread.fan(--i).m_state = HostThread::ThreadActive ;
  }
}

/** \brief  This thread waits for each fan-in thread to deactivate.  */
void HostThreadWorker::fanin_deactivation( HostThread & thread ) const
{
  // The 'wait' function repeatedly polls the 'thread' state
  // which may reside in a different NUMA region.
  // Thus the fan is intra-node followed by inter-node
  // to minimize inter-node memory access.

  for ( unsigned i = 0 ; i < thread.fan_count() ; ++i ) {
    host_wait( & thread.fan(i).m_state , HostThread::ThreadActive );
  }
}

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
    host_wait( & thread.m_state , HostThread::ThreadActive );
  }

  // Execute on the master thread,
  worker.execute_on_thread( m_master_thread );

  // Worker threads are returned to the ThreadInactive state.
  m_worker = NULL ;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

HostInternal::~HostInternal()
{
  HostThread::clear_thread(0);
}

HostInternal::HostInternal()
  : m_worker_block()

  , m_gang_capacity( 1 )
  , m_worker_capacity( 0 )

  , m_cache_line_size( 64 /* default guess */ )
  , m_thread_count( 1 )
  , m_gang_count( 1 )
  , m_worker_count( 1 )
  , m_reduce_scratch_size( 0 )
{
  m_worker = NULL ;

  if ( ! is_master_thread() ) {
    KokkosArray::Impl::throw_runtime_exception( std::string("KokkosArray::Impl::HostInternal FAILED : not initialized on the master thread") );
  }

  // Initialize thread pool with a single master thread:
  HostThread::set_thread( 0 , & m_master_thread );
  HostThread::set_thread_relationships();
}

void HostInternal::print_configuration( std::ostream & s ) const
{
  s << "KokkosArray::Host thread capacity( "
    << m_gang_capacity << " x " << m_worker_capacity
    << " ) used( "
    << m_gang_count << " x " << m_worker_count
    << " )"
    << std::endl ;
}

//----------------------------------------------------------------------------

inline
void HostInternal::resize_reduce_thread( HostThread & thread ) const
{
  if ( thread.m_reduce ) {
    host_decrement_not_thread_safe( thread.m_reduce );
    thread.m_reduce = 0 ;
  }

  if ( m_reduce_scratch_size ) {

    thread.m_reduce = host_allocate_not_thread_safe( "reduce_scratch_space" , typeid(unsigned char) , 1 , m_reduce_scratch_size );

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
void HostInternal::resize_reduce_scratch( unsigned size )
{
  if ( 0 == size || m_reduce_scratch_size < size ) {
    // Round up to cache line size:
    const unsigned rem = size % m_cache_line_size ;

    if ( rem ) size += m_cache_line_size - rem ;
  
    m_reduce_scratch_size = size ;

    const HostWorkerResizeReduce work ;

    execute_serial( work );
  }
}

void * HostInternal::reduce_scratch() const
{
  return m_master_thread.reduce_data();
}

void host_resize_scratch_reduce( unsigned size )
{ HostInternal::singleton().resize_reduce_scratch( size ); }

void * host_scratch_reduce()
{ return HostInternal::singleton().reduce_scratch(); }

//----------------------------------------------------------------------------

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

    HostThread * const thread = HostThread::get_thread( m_thread_count );

    if ( 0 != thread ) {

      m_master_thread.m_state = HostThread::ThreadInactive ;

      thread->m_state = HostThread::ThreadTerminating ;

      host_wait( & m_master_thread.m_state , HostThread::ThreadInactive );

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

  if ( bind_thread( thread_rank ) ) {

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
        host_wait( & this_thread.m_state , HostThread::ThreadInactive );
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

  bool ok_spawn_threads = true ;

  for ( unsigned rank = gang_count * worker_count ; ok_spawn_threads && 0 < --rank ; ) {

    m_master_thread.m_state = HostThread::ThreadInactive ;

    // Spawn thread executing the 'driver' function.
    ok_spawn_threads = spawn( rank );

    if ( ok_spawn_threads ) {

      // Thread spawned, wait for thread to activate:
      host_wait( & m_master_thread.m_state , HostThread::ThreadInactive );

      // Check if the thread initialized and bound correctly:
      ok_spawn_threads = HostThread::ThreadActive == m_master_thread.m_state ;

      if ( ok_spawn_threads ) { // Wait for spawned thread to deactivate
        host_wait( & HostThread::get_thread(rank)->m_state , HostThread::ThreadActive );
      }
    }
  }

  // Bind the process thread as thread_rank == 0
  if ( ok_spawn_threads ) { ok_spawn_threads = bind_thread( 0 ); }

  m_worker = NULL ;

  return ok_spawn_threads ;
}

void HostInternal::initialize( const unsigned gang_count ,
                               const unsigned worker_count )
{
  const bool ok_inactive     = 1 == m_thread_count ;
  const bool ok_gang_count   = gang_count <= m_gang_capacity ;
  const bool ok_worker_count = ( 0 == m_worker_capacity ||
                                 worker_count <= m_worker_capacity );

  // Only try to spawn threads if input is valid.

  const bool ok_spawn_threads =
    ( ok_inactive &&
      ok_gang_count && ok_worker_count &&
      1 < gang_count * worker_count )
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
          << m_gang_capacity << " x " << m_worker_capacity
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
                       const Host::size_type gang_worker_count )
{
  Impl::HostInternal::singleton().initialize( gang_count , gang_worker_count );
}

void Host::print_configuration( std::ostream & s )
{ Impl::HostInternal::singleton().print_configuration(s); }

Host::size_type Host::detect_gang_capacity()
{ return Impl::HostInternal::singleton().m_gang_capacity ; }

Host::size_type Host::detect_gang_worker_capacity()
{ return Impl::HostInternal::singleton().m_worker_capacity ; }

Host::size_type Host::detect_cache_line_size()
{ return Impl::HostInternal::singleton().m_cache_line_size ; }

/*--------------------------------------------------------------------------*/

} // namespace KokkosArray


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
#include <KokkosArray_hwloc.hpp>
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
namespace {

class HostWorkerBlock : public HostThreadWorker {
public:

  void execute_on_thread( HostThread & thread ) const 
  {
    host_thread_lock();
    host_thread_unlock();
    end_barrier( thread );
  }

  HostWorkerBlock()  {}
  ~HostWorkerBlock() {}
};

const HostWorkerBlock worker_block ;

std::pair<unsigned,unsigned> host_thread_coordinates[ HostThread::max_thread_count ] ;

//----------------------------------------------------------------------------

inline
void thread_mapping( const unsigned thread_rank ,
                     const unsigned gang_count,
                     const unsigned worker_count ,
                     std::pair<unsigned,unsigned> & coordinate )
{
  const unsigned gang_rank   = thread_rank / worker_count ;
  const unsigned worker_rank = thread_rank % worker_count ;

  const std::pair<unsigned,unsigned> core_topo = hwloc::get_core_topology();

  { // Distribute gangs amont NUMA regions:
    // gang_count = k * bin + ( #NUMA - k ) * ( bin + 1 )
    const unsigned bin  = gang_count / core_topo.first ;
    const unsigned bin1 = bin + 1 ;
    const unsigned k    = core_topo.first * bin1 - gang_count ;
    const unsigned part = k * bin ;

    coordinate.first = ( gang_rank < part )
                     ? ( gang_rank / bin )
                     : ( k + ( gang_rank - part ) / bin1 );
  }

  { // Distribute workers to cores:
    // worker_count = k * bin + ( (#CORE/NUMA) - k ) * ( bin + 1 )
    const unsigned bin  = worker_count / core_topo.second ;
    const unsigned bin1 = bin + 1 ;
    const unsigned k    = core_topo.second * bin1 - worker_count ;
    const unsigned part = k * bin ;

    coordinate.second = ( worker_rank < part )
                      ? ( worker_rank / bin )
                      : ( k + ( worker_rank - part ) / bin1 );
  }
}

} // namespace

//----------------------------------------------------------------------------

void HostInternal::activate_threads()
{
  if ( 1 < m_thread_count ) {
    // Activate threads to call 'm_worker.execute_on_thread'
    for ( unsigned i = m_thread_count ; 0 < i ; ) {
      HostThread::get_thread(--i)->m_state = HostThread::ThreadActive ;
    }
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

  // Wait for threads to complete:
  worker.end_barrier( m_master_thread );

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

  for ( unsigned rank = m_thread_count ; 1 < rank ; ) {
    --rank ;

    for ( unsigned j = 0 ; j < m_thread_count ; ++j ) {

      HostThread & thread = * HostThread::get_thread(j);

      if ( thread.rank() == rank ) {
        thread.m_state = HostThread::ThreadActive ;
        host_thread_wait( & thread.m_state , HostThread::ThreadActive );
        break ;
      }
    }
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
{
  m_worker = NULL ;

  if ( ! is_master_thread() ) {
    KokkosArray::Impl::throw_runtime_exception( std::string("KokkosArray::Impl::HostInternal FAILED : not initialized on the master thread") );
  }
}

HostInternal & HostInternal::singleton()
{
  static HostInternal self ; return self ;
}

void HostInternal::print_configuration( std::ostream & s ) const
{
  const std::pair<unsigned,unsigned> core_topo = hwloc::get_core_topology();
  const unsigned core_size = hwloc::get_core_capacity();

  s << "NUMA[" << core_topo.first << "]"
    << " CORE[" << core_topo.second << "]"
    << " PU[" << core_size << "]"
    << std::endl ;

  for ( unsigned r = 0 ; r < m_thread_count ; ++r ) {
    HostThread * const thread = HostThread::get_thread( r );

    std::pair<unsigned,unsigned> coord ;

    thread_mapping( r , thread->gang_count() ,
                        thread->worker_count() ,
                        coord );

    s << "  thread[ " << thread->rank() << " : ("
      << thread->gang_rank()
      << "," << thread->worker_rank() << ") ] bound to core("
      << coord.first << "," << coord.second << ")"
      << std::endl ;
  }
}

//----------------------------------------------------------------------------

inline
void * HostInternal::reduce_scratch() const
{
  return m_master_thread.reduce_data();
}

//----------------------------------------------------------------------------

void HostInternal::finalize()
{
  verify_inactive("finalize()");

  Host::resize_reduce_scratch(0);

  // Release and clear worker threads:

  HostThread::clear_thread_relationships();

  for ( unsigned r = 0 ; r < HostThread::max_thread_count ; ++r ) {

    HostThread * thread = HostThread::get_thread( r );

    if ( & m_master_thread != thread && 0 != thread ) {

      m_master_thread.m_state = HostThread::ThreadInactive ;

      thread->m_state = HostThread::ThreadTerminating ;

      thread = 0 ; // The '*thread' object is now invalid

      // Wait for '*thread' to terminate:
      host_thread_wait( & m_master_thread.m_state , HostThread::ThreadInactive );
    }

    HostThread::clear_thread(r);
  }

  m_thread_count = 1 ;

  hwloc::unbind_this_thread();
}

//----------------------------------------------------------------------------
// Bind this thread to one of the requested cores and remove that request.

unsigned HostInternal::bind_host_thread()
{
  const std::pair<unsigned,unsigned> current = hwloc::get_thread_coordinate();

  unsigned i = m_thread_count ;

  if ( m_thread_count ) {

    // Match one of the requests:
    for ( i = 0 ; i < m_thread_count && current != host_thread_coordinates[i] ; ++i );

    if ( m_thread_count == i ) {
      // Match the NUMA request:
      for ( i = 0 ; i < m_thread_count && current.first != host_thread_coordinates[i].first ; ++i );
    }

    if ( m_thread_count == i ) {
      // Match any unclaimed request:
      for ( i = 0 ; i < m_thread_count && ~0u == host_thread_coordinates[i].first  ; ++i );
    }

    if ( i < m_thread_count ) {
      if ( ! hwloc::bind_this_thread( host_thread_coordinates[i] ) ) i = m_thread_count ;
    }

    if ( i < m_thread_count ) {
      host_thread_coordinates[i].first  = ~0u ;
      host_thread_coordinates[i].second = ~0u ;
    }
  }

  return i ;
}

// Driver for each created thread

void HostInternal::driver()
{
  HostInternal & host = Impl::HostInternal::singleton();

  // Bind this thread to a unique processing unit.

  const unsigned thread_entry = host.bind_host_thread();

  if ( thread_entry < host.m_thread_count ) {

    HostThread this_thread ;

    this_thread.m_state = HostThread::ThreadActive ;

    HostThread::set_thread( thread_entry , & this_thread );

    // Inform master thread that spawning and binding succeeded.
    host.m_master_thread.m_state = HostThread::ThreadActive ;

    // Work loop:

    while ( HostThread::ThreadActive == this_thread.m_state ) {

      // Perform the work:
      host.m_worker->execute_on_thread( this_thread );

      // Deactivate this thread:
      this_thread.m_state = HostThread::ThreadInactive ;

      // Wait to be activated or terminated:
      host_thread_wait( & this_thread.m_state , HostThread::ThreadInactive );
    }
  }

  // Notify master thread that this thread has terminated.

  host.m_master_thread.m_state = HostThread::ThreadTerminating ;
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

  const unsigned thread_count = gang_count * worker_count ;

  // Reserve master thread's coordinates:
  std::pair<unsigned,unsigned> master_core = hwloc::get_thread_coordinate();

  // Define coordinates for pinning threads:

  for ( unsigned r = 0 ; r < thread_count ; ++r ) {
    thread_mapping( r , gang_count , worker_count , host_thread_coordinates[r] );
  }

  {
    unsigned i = 0 ;
    for ( ; i < m_thread_count && master_core != host_thread_coordinates[i] ; ++i );

    if ( i == m_thread_count ) {
      for ( i = 0 ; i < m_thread_count &&
                    master_core.first != host_thread_coordinates[i].first ; ++i );
    }

    if ( i == m_thread_count ) i = 0 ;

    host_thread_coordinates[i] = host_thread_coordinates[0] ;
    host_thread_coordinates[0] = std::pair<unsigned,unsigned>(~0u,~0u);
  }

  bool ok_spawn_threads = true ;

  m_worker = & worker_block ;

  for ( unsigned i = 1 ; i < thread_count && ok_spawn_threads ; ++i ) {

    m_master_thread.m_state = HostThread::ThreadInactive ;

    // Spawn thread executing the 'driver' function.
    ok_spawn_threads = spawn();

    if ( ok_spawn_threads ) {

      // Thread spawned, wait for thread to activate:
      host_thread_wait( & m_master_thread.m_state , HostThread::ThreadInactive );

      // Check if the thread initialized and bound correctly:
      ok_spawn_threads = HostThread::ThreadActive == m_master_thread.m_state ;
    }
  }

  // All successfully spawned threads are in the list...

  unsigned count_verify = 0 ;

  for ( unsigned entry = 0 ; entry < thread_count ; ++entry ) {
    HostThread * const thread = HostThread::get_thread(entry);

    if ( thread ) {
      host_thread_wait( & thread->m_state , HostThread::ThreadActive );
      ++count_verify ;
    }
  }

  ok_spawn_threads = thread_count == 1 + count_verify ;

  if ( ok_spawn_threads ) {
    hwloc::bind_this_thread( master_core );
    HostThread::set_thread( 0 , & m_master_thread );
  }

  m_worker = NULL ;

  return ok_spawn_threads ;
}

void HostInternal::initialize( const unsigned gang_count ,
                               const unsigned worker_count )
{
  const bool ok_inactive = 1 == m_thread_count &&
                           0 == HostThread::get_thread_count();
  bool ok_spawn_threads = true ;

  // Only try to spawn threads if input is valid.

  if ( ok_inactive && 1 < gang_count * worker_count ) {

    if ( ok_spawn_threads ) {

      m_thread_count = gang_count * worker_count ;

      ok_spawn_threads = spawn_threads( gang_count , worker_count );

      if ( ok_spawn_threads ) {
        // All threads spawned, initialize the master-thread fan-in

        for ( unsigned g = 0 , r = 0 ; g < gang_count ; ++g ) {
        for ( unsigned w = 0 ;         w < worker_count ; ++w , ++r ) {

          HostThread & thread = * HostThread::get_thread(r);

          thread.m_gang_tag    = g ;
          thread.m_worker_tag  = w ;
        }}

        if ( 0 != m_master_thread.m_gang_tag ) {
          ok_spawn_threads = false ;
        }
        else if ( 0 != m_master_thread.m_worker_tag ) {
          // The master thread must be (0,0)
          // Swap whatever thread happened to land there.
          HostThread & thread = * HostThread::get_thread(0);

          thread.m_worker_tag = m_master_thread.m_worker_tag ;
          m_master_thread.m_worker_tag = 0 ;
        }

        HostThread::set_thread_relationships();
      }
      else {
        finalize();
      }
    }
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

  // Initial reduction allocation:
  Host::resize_reduce_scratch( 4096 );
}

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

//----------------------------------------------------------------------------

namespace Impl {
namespace {

struct HostThreadResizeReduce : public HostThreadWorker {

  const unsigned reduce_size ;

  HostThreadResizeReduce( unsigned size )
    : reduce_size(size)
    { HostThreadWorker::execute_serial(); }

  void execute_on_thread( HostThread & thread ) const
    {
      thread.resize_reduce( reduce_size );
      end_barrier( thread );
    }
};

}
} // namespace Impl

void Host::resize_reduce_scratch( unsigned size )
{
  static unsigned m_reduce_size = 0 ;

  const unsigned rem = size % HostSpace::MEMORY_ALIGNMENT ;

  if ( rem ) size += HostSpace::MEMORY_ALIGNMENT - rem ;

  if ( ( 0 == size ) || ( m_reduce_size < size ) ) {

    Impl::HostThreadResizeReduce tmp(size);

    m_reduce_size = size ;
  }
}

void * Host::root_reduce_scratch()
{ return Impl::HostInternal::singleton().reduce_scratch(); }

//----------------------------------------------------------------------------

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
  const std::pair<unsigned,unsigned> core_topo = hwloc::get_core_topology();

  return core_topo.first ;
}

Host::size_type Host::detect_gang_worker_capacity()
{
  const std::pair<unsigned,unsigned> core_topo = hwloc::get_core_topology();
  const unsigned core_size = hwloc::get_core_capacity();

  return core_topo.second * core_size ;
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

    Impl::worker_block.end_barrier( h.m_master_thread );

    h.m_worker = NULL ;

    is_ready = true ;
  }

  return is_ready ;
}

/*--------------------------------------------------------------------------*/

} // namespace KokkosArray


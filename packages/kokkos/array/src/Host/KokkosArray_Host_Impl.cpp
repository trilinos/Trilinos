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

void host_thread_mapping( const std::pair<unsigned,unsigned> gang_topo ,
                          const std::pair<unsigned,unsigned> core_use ,
                          const std::pair<unsigned,unsigned> core_topo ,
                          const std::pair<unsigned,unsigned> master_coord ,
                                std::pair<unsigned,unsigned> thread_coord[] )
{
  const unsigned thread_count = gang_topo.first * gang_topo.second ;
  const unsigned core_base    = core_topo.second - core_use.second ;

  for ( unsigned thread_rank = 0 , gang_rank = 0 ; gang_rank < gang_topo.first ; ++gang_rank ) {
  for ( unsigned worker_rank = 0 ; worker_rank < gang_topo.second ; ++worker_rank , ++thread_rank ) {

    unsigned gang_in_numa_count = 0 ;
    unsigned gang_in_numa_rank  = 0 ;

    { // Distribute gangs among NUMA regions:
      // gang_count = k * bin + ( #NUMA - k ) * ( bin + 1 )
      const unsigned bin  = gang_topo.first / core_use.first ;
      const unsigned bin1 = bin + 1 ;
      const unsigned k    = core_use.first * bin1 - gang_topo.first ;
      const unsigned part = k * bin ;

      if ( gang_rank < part ) {
        thread_coord[ thread_rank ].first = gang_rank / bin ;
        gang_in_numa_rank  = gang_rank % bin ;
        gang_in_numa_count = bin ;
      }
      else {
        thread_coord[ thread_rank ].first = k + ( gang_rank - part ) / bin1 ;
        gang_in_numa_rank  = ( gang_rank - part ) % bin1 ;
        gang_in_numa_count = bin1 ;
      }
    }

    { // Distribute workers to cores within this NUMA region:
      // worker_in_numa_count = k * bin + ( (#CORE/NUMA) - k ) * ( bin + 1 )
      const unsigned worker_in_numa_count = gang_in_numa_count * gang_topo.second ;
      const unsigned worker_in_numa_rank  = gang_in_numa_rank  * gang_topo.second + worker_rank ;

      const unsigned bin  = worker_in_numa_count / core_use.second ;
      const unsigned bin1 = bin + 1 ;
      const unsigned k    = core_use.second * bin1 - worker_in_numa_count ;
      const unsigned part = k * bin ;

      thread_coord[ thread_rank ].second = core_base +
        ( ( worker_in_numa_rank < part )
          ? ( worker_in_numa_rank / bin )
          : ( k + ( worker_in_numa_rank - part ) / bin1 ) );
    }
  }}

  // The master core should be thread #0 so rotate all coordinates accordingly ...

  const std::pair<unsigned,unsigned> offset
    ( ( thread_coord[0].first  < master_coord.first  ? master_coord.first  - thread_coord[0].first  : 0 ) ,
      ( thread_coord[0].second < master_coord.second ? master_coord.second - thread_coord[0].second : 0 ) );

  for ( unsigned i = 0 ; i < thread_count ; ++i ) {
    thread_coord[i].first  = ( thread_coord[i].first + offset.first ) % core_use.first ;
    thread_coord[i].second = core_base + ( thread_coord[i].second + offset.second - core_base ) % core_use.second ;
  }

#if 0

  std::cout << "KokkosArray::Host thread_mapping" << std::endl ;

  for ( unsigned g = 0 , t = 0 ; g < gang_topo.first ; ++g ) {
    std::cout << "  gang[" << g 
              << "] on numa[" << thread_coord[t].first
              << "] cores(" ;
    for ( unsigned w = 0 ; w < gang_topo.second ; ++w , ++t ) {
      std::cout << " " << thread_coord[t].second ;
    }
    std::cout << " )" << std::endl ;
  }

#endif

}

} // namespace Impl
} // namespace KokkosArray

/*--------------------------------------------------------------------------*/

namespace KokkosArray {
namespace {

class HostWorkerBlock : public Impl::HostThreadWorker {
public:

  void execute_on_thread( Impl::HostThread & thread ) const 
  {
    Impl::host_thread_lock();
    Impl::host_thread_unlock();
    end_barrier( thread );
  }

  HostWorkerBlock()  {}
  ~HostWorkerBlock() {}
};

//----------------------------------------------------------------------------

const HostWorkerBlock s_worker_block ;
Impl::HostThread      s_master_thread ;
unsigned              s_host_thread_count = 0 ;

std::pair<unsigned,unsigned>
  s_host_thread_coord[ Impl::HostThread::max_thread_count ] ;

const Impl::HostThreadWorker * volatile s_current_worker = 0 ;

//----------------------------------------------------------------------------

struct Sentinel {
  Sentinel() {}
  ~Sentinel();
};

Sentinel::~Sentinel()
{
  if ( 1 < s_master_thread.count() ) {
    std::cerr << "KokkosArray::Host ERROR : initialized but not finalized"
              << std::endl ;
  }
}

//----------------------------------------------------------------------------
// Bind this thread to one of the requested cores and remove that request.

unsigned bind_host_thread()
{
  const std::pair<unsigned,unsigned> current = hwloc::get_this_thread_coordinate();

  unsigned i = 0 ;

  // Match one of the requests:
  for ( i = 0 ; i < s_host_thread_count && current != s_host_thread_coord[i] ; ++i );

  if ( s_host_thread_count == i ) {
    // Match the NUMA request:
    for ( i = 0 ; i < s_host_thread_count && current.first != s_host_thread_coord[i].first ; ++i );
  }

  if ( s_host_thread_count == i ) {
    // Match any unclaimed request:
    for ( i = 0 ; i < s_host_thread_count && ~0u == s_host_thread_coord[i].first  ; ++i );
  }

  if ( i < s_host_thread_count ) {
    if ( ! hwloc::bind_this_thread( s_host_thread_coord[i] ) ) i = s_host_thread_count ;
  }

  if ( i < s_host_thread_count ) {

#if 0

    if ( current != s_host_thread_coord[i] ) {
      std::cout << "  rebinding thread[" << i << "] from ("
                << current.first << "," << current.second
                << ") to ("
                << s_host_thread_coord[i].first << ","
                << s_host_thread_coord[i].second
                << ")" << std::endl ;
    }

#endif

    s_host_thread_coord[i].first  = ~0u ;
    s_host_thread_coord[i].second = ~0u ;
  }

  return i ;
}


bool spawn_threads( const std::pair<unsigned,unsigned> gang_topo ,
                          std::pair<unsigned,unsigned> core_use )
{
  // If the process is bound to a particular node
  // then only use cores belonging to that node.
  // Otherwise use all nodes and all their cores.

  // Define coordinates for pinning threads:
  const unsigned thread_count = gang_topo.first * gang_topo.second ;

  bool ok_spawn_threads = true ;

  if ( 1 < thread_count ) {
    const std::pair<unsigned,unsigned> core_topo    = hwloc::get_core_topology();
    const std::pair<unsigned,unsigned> master_coord = hwloc::get_this_thread_coordinate();

    if ( 0 == core_use.first  || core_topo.first  < core_use.first ||
         0 == core_use.second || core_topo.second < core_use.second ) {
      core_use = core_topo ;
    }

    Impl::host_thread_mapping( gang_topo , core_use , core_topo , master_coord , s_host_thread_coord );

    // Reserve the master thread coordinate for the mater thread:

    const std::pair<unsigned,unsigned> master_core = s_host_thread_coord[0] ;

    s_host_thread_coord[0] = std::pair<unsigned,unsigned>(~0u,~0u);

    s_current_worker = & s_worker_block ;
    s_host_thread_count = thread_count ;

    for ( unsigned i = 1 ; i < thread_count && ok_spawn_threads ; ++i ) {

      s_master_thread.m_state = Impl::HostThread::ThreadInactive ;

      // Spawn thread executing the 'driver' function.
      ok_spawn_threads = Impl::host_thread_spawn();

      // Thread spawned.  Thread will inform me of its condition as:
      //   Active     = success and has entered its 'HostThread' data in the lookup table.
      //   Terminate  = failure

      if ( ok_spawn_threads ) {
        ok_spawn_threads =
          Impl::HostThread::ThreadActive ==
            host_thread_wait( & s_master_thread.m_state , Impl::HostThread::ThreadInactive );
      }
    }

    // All successfully spawned threads are in the list as [1,thread_count)
    // Wait for them to deactivate before zeroing the worker

    for ( unsigned rank = 1 ; rank < thread_count ; ++rank ) {
      Impl::HostThread * const th = Impl::HostThread::get_thread(rank);
      if ( th ) {
        host_thread_wait( & th->m_state , Impl::HostThread::ThreadActive );
      }
    }

    s_current_worker = NULL ;
    s_host_thread_count = 0 ;

    if ( ok_spawn_threads ) {

      hwloc::bind_this_thread( master_core );

      Impl::HostThread::set_thread( 0 , & s_master_thread );

      // All threads spawned, set the fan-in relationships

      for ( unsigned g = 0 , r = 0 ; g < gang_topo.first ; ++g ) {
      for ( unsigned w = 0 ;         w < gang_topo.second ; ++w , ++r ) {
        Impl::HostThread::get_thread(r)->set_topology( r , thread_count ,
                                                       g , gang_topo.first ,
                                                       w , gang_topo.second );
      }}

      Impl::HostThread::set_thread_relationships();
    }
  }

  return ok_spawn_threads ;
}

//----------------------------------------------------------------------------

inline
void activate_threads()
{
  for ( unsigned i = s_master_thread.count() ; 0 < --i ; ) {
    Impl::HostThread::get_thread(i)->m_state = Impl::HostThread::ThreadActive ;
  }
}

void verify_inactive( const char * const method )
{
  if ( NULL != s_current_worker ) {

    std::ostringstream msg ;
    msg << method << " FAILED: " ;

    if ( & s_worker_block == s_current_worker ) {
      msg << "Host is blocked" ;
    }
    else {
      msg << "Functor is executing." ;
    }
    KokkosArray::Impl::throw_runtime_exception( msg.str() );
  }
}

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

// Driver called by each created thread

void host_thread_driver()
{
  // Bind this thread to a unique processing unit.

  const unsigned thread_rank = bind_host_thread();

  if ( thread_rank < s_host_thread_count ) {

    HostThread this_thread ;

    this_thread.m_state = HostThread::ThreadActive ;

    HostThread::set_thread( thread_rank , & this_thread );

    // Inform master thread that spawning and binding succeeded.
    s_master_thread.m_state = HostThread::ThreadActive ;

    // Work loop:

    while ( HostThread::ThreadActive == this_thread.m_state ) {

      // Perform the work:
      s_current_worker->execute_on_thread( this_thread );

      // Deactivate this thread:
      this_thread.m_state = HostThread::ThreadInactive ;

      // Wait to be activated or terminated:
      host_thread_wait( & this_thread.m_state , HostThread::ThreadInactive );
    }
  }

  // Notify master thread that this thread has terminated.

  s_master_thread.m_state = HostThread::ThreadTerminating ;
}

//----------------------------------------------------------------------------

void HostThreadWorker::execute() const
{
  verify_inactive("KokkosArray::Impl::HostThreadWorker::execute()");

  // Worker threads are verified to be in the ThreadInactive state.

  s_current_worker = this ;

  activate_threads();

  // Execute on the master thread,
  execute_on_thread( s_master_thread );

  // Wait for threads to complete:
  end_barrier( s_master_thread );

  // Worker threads are returned to the ThreadInactive state.
  s_current_worker = NULL ;
}

void HostThreadWorker::execute_serial() const
{
  verify_inactive("KokkosArray::Impl::HostThreadWorker::execute_serial()");

  // Worker threads are verified to be in the ThreadInactive state.

  s_current_worker = this ;

  for ( unsigned rank = s_master_thread.count() ; 1 < rank ; ) {
    --rank ;

    HostThread & thread = * HostThread::get_thread(rank);

    thread.m_state = HostThread::ThreadActive ;

    host_thread_wait( & thread.m_state , HostThread::ThreadActive );
  }

  // Execute on the master thread last.
  execute_on_thread( s_master_thread );

  // Worker threads are returned to the ThreadInactive state.
  s_current_worker = NULL ;
}

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace {

struct HostThreadResizeReduce : public Impl::HostThreadWorker {

  const unsigned reduce_size ;

  HostThreadResizeReduce( unsigned size )
    : reduce_size(size)
    { HostThreadWorker::execute_serial(); }

  void execute_on_thread( Impl::HostThread & thread ) const
    {
      thread.resize_reduce( reduce_size );
      end_barrier( thread );
    }
};

}

void Host::resize_reduce_scratch( unsigned size )
{
  static unsigned m_reduce_size = 0 ;

  const unsigned rem = size % HostSpace::MEMORY_ALIGNMENT ;

  if ( rem ) size += HostSpace::MEMORY_ALIGNMENT - rem ;

  if ( ( 0 == size ) || ( m_reduce_size < size ) ) {

    HostThreadResizeReduce tmp(size);

    m_reduce_size = size ;
  }
}

void * Host::root_reduce_scratch()
{ return s_master_thread.reduce_data(); }

//----------------------------------------------------------------------------

void Host::finalize()
{
  verify_inactive("KokkosArray::Host::finalize()");

  Host::resize_reduce_scratch(0);

  // Release and clear worker threads:

  Impl::HostThread::clear_thread_relationships();

  for ( unsigned r = 0 ; r < Impl::HostThread::max_thread_count ; ++r ) {

    Impl::HostThread * thread = Impl::HostThread::get_thread( r );

    if ( & s_master_thread != thread && 0 != thread ) {

      s_master_thread.m_state = Impl::HostThread::ThreadInactive ;

      thread->m_state = Impl::HostThread::ThreadTerminating ;

      thread = 0 ; // The '*thread' object is now invalid

      // Wait for '*thread' to terminate:
      host_thread_wait( & s_master_thread.m_state , Impl::HostThread::ThreadInactive );
    }

    Impl::HostThread::clear_thread(r);
  }

  hwloc::unbind_this_thread();
}

void Host::initialize( const std::pair<unsigned,unsigned> gang_topo ,
                       const std::pair<unsigned,unsigned> core_topo )
{
  static const Sentinel sentinel ;

  const bool ok_inactive = Impl::host_thread_is_master() &&
                           1 == s_master_thread.count() &&
                           0 == Impl::HostThread::get_thread(0);
  bool ok_spawn_threads = true ;

  // Only try to spawn threads if input is valid.

  if ( ok_inactive && 1 < gang_topo.first * gang_topo.second ) {

    ok_spawn_threads = spawn_threads( gang_topo , core_topo );

    if ( ! ok_spawn_threads ) {
      finalize();
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

void Host::print_configuration( std::ostream & s )
{
  const std::pair<unsigned,unsigned> core_topo = hwloc::get_core_topology();
  const unsigned core_size = hwloc::get_core_capacity();

  s << "hwloc { NUMA[" << core_topo.first << "]"
    << " CORE[" << core_topo.second << "]"
    << " PU[" << core_size << "] } "
    << "threadpool { GANG[" << s_master_thread.gang_count() << "]"
    << " WORKER[" << s_master_thread.worker_count() << "] }"
    << std::endl ;
}

//----------------------------------------------------------------------------

bool Host::sleep()
{
  const bool is_ready   = NULL == s_current_worker ;
        bool is_blocked = & s_worker_block == s_current_worker ;

  if ( is_ready ) {
    Impl::host_thread_lock();

    s_current_worker = & s_worker_block ;

    // Activate threads so that they will proceed from
    // spinning state to being blocked on the mutex.

    activate_threads();

    is_blocked = true ;
  }

  return is_blocked ;
}

bool Host::wake()
{
  const bool is_blocked = & s_worker_block == s_current_worker ;
        bool is_ready   = NULL == s_current_worker ;

  if ( is_blocked ) {
    Impl::host_thread_unlock();

    s_worker_block.end_barrier( s_master_thread );

    s_current_worker = NULL ;

    is_ready = true ;
  }

  return is_ready ;
}

/*--------------------------------------------------------------------------*/

} // namespace KokkosArray


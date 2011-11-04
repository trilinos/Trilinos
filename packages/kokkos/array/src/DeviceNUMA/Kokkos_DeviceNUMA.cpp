/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2011 Sandia Corporation
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

/*--------------------------------------------------------------------------*/
/* Kokkos interfaces */

#include <Kokkos_DeviceNUMA.hpp>
#include <impl/Kokkos_MemoryInfo.hpp>

/*--------------------------------------------------------------------------*/
/* Standard 'C' libraries */
#include <stdlib.h>

/* Standard 'C++' libraries */
#include <iostream>
#include <stdexcept>
#include <sstream>

/*--------------------------------------------------------------------------*/
/* Third Party Libraries */

/* Hardware locality library: http://www.open-mpi.org/projects/hwloc/ */
#include <hwloc.h>

#define  REQURED_HWLOC_API_VERSION  0x000010300

#if HWLOC_API_VERSION < REQUIRED_HWLOC_VERSION
#error "Requires  http://www.open-mpi.org/projects/hwloc/  Version 1.3 or greater"
#endif

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

void device_numa_thread_lock();
void device_numa_thread_unlock();
bool device_numa_thread_spawn( DeviceNUMAThread * );

//----------------------------------------------------------------------------

class DeviceNUMAWorkerBlock : public DeviceNUMAWorker {
public:
  void execute_on_thread( DeviceNUMAThread & ) const ;

  DeviceNUMAWorkerBlock()  {}
  ~DeviceNUMAWorkerBlock() {}
};

void DeviceNUMAWorkerBlock::execute_on_thread( DeviceNUMAThread & this_thread ) const
{
  device_numa_thread_lock();
  device_numa_thread_unlock();

  this_thread.barrier();
}

//----------------------------------------------------------------------------

class DeviceNUMAInternal {
public:

  enum { THREAD_COUNT_MAX = 255 };

  typedef DeviceNUMA::size_type size_type ;

  hwloc_topology_t      m_host_topology ;
  DeviceNUMAWorkerBlock m_worker_block ;
  size_type             m_node_count ;
  size_type             m_cpu_per_node ;
  size_type             m_thread_count ;
  size_type             m_thread_per_node ;
  DeviceNUMAThread      m_thread[ THREAD_COUNT_MAX + 1 ];

  const DeviceNUMAWorker * volatile m_worker ;

  ~DeviceNUMAInternal();
  DeviceNUMAInternal();

  void verify_inactive( const char * const method ) const ;

  void initialize( DeviceNUMA::UtilizationOptions );
  void finalize();
  void block();
  void unblock();

  bool cpubind( size_type node_rank ) const ;

  void execute( const DeviceNUMAWorker & worker );

  static DeviceNUMAInternal & singleton();
};

DeviceNUMAInternal & DeviceNUMAInternal::singleton()
{
  static DeviceNUMAInternal self ;
  return self ;
}

DeviceNUMAInternal::~DeviceNUMAInternal()
{
  hwloc_topology_destroy( m_host_topology );
}

DeviceNUMAInternal::DeviceNUMAInternal()
{
  m_node_count      = 0 ;
  m_cpu_per_node    = 0 ;
  m_thread_count    = 0 ;
  m_thread_per_node = 0 ;
  m_worker          = NULL ;

  hwloc_topology_init( & m_host_topology );
  hwloc_topology_load( m_host_topology );

  const int hwloc_depth =
    hwloc_get_type_depth( m_host_topology , HWLOC_OBJ_NODE );

  if ( HWLOC_TYPE_DEPTH_UNKNOWN != hwloc_depth ) {

    m_node_count =
      hwloc_get_nbobjs_by_depth( m_host_topology , hwloc_depth );

    bool ok = 0 < m_node_count ;

    for ( size_type i = 0 ; i < m_node_count && ok ; ++i ) {

      const hwloc_obj_t node =
        hwloc_get_obj_by_type( m_host_topology , HWLOC_OBJ_NODE , i );

      const size_type count = hwloc_bitmap_weight( node->allowed_cpuset );

      if ( 0 == m_cpu_per_node ) { m_cpu_per_node = count ; }

      ok = count == m_cpu_per_node ;
    }

    if ( ok ) {
      ok = cpubind( 0 ); // Bind this master thread
    }

    if ( ! ok ) {
      m_cpu_per_node = 0 ;
      m_node_count = 0 ;
    }
  }
}

// Bind the current thread to the location for the input thread rank

bool DeviceNUMAInternal::cpubind( DeviceNUMA::size_type node_rank ) const
{
  const hwloc_obj_t node =
    hwloc_get_obj_by_type( m_host_topology, HWLOC_OBJ_NODE, node_rank );

  return -1 != hwloc_set_cpubind( m_host_topology , node->allowed_cpuset ,
                                  HWLOC_CPUBIND_THREAD | HWLOC_CPUBIND_STRICT );
}

//----------------------------------------------------------------------------
// Driver for each created thread

void DeviceNUMAThread::driver()
{
  const DeviceNUMAInternal & pool = DeviceNUMAInternal::singleton();

  //------------------------------------
  // this_thread is in the Active state.
  // The master thread is waiting for this_thread to leave the Active state.
  //------------------------------------
  // Migrate myself to the proper NUMA node and run

  if ( pool.cpubind( rank() / pool.m_thread_per_node ) ) {

    while ( ThreadActive == m_state ) {

      // When the work is complete the state will be Inactive or Terminate
      pool.m_worker->execute_on_thread( *this );

      // If this_thread is in the Inactive state then wait for activation.
      wait( ThreadInactive );
    }
  }

  set( ThreadNull );
}

//----------------------------------------------------------------------------

void DeviceNUMAInternal::verify_inactive( const char * const method ) const
{
  if ( NULL != m_worker ) {

    std::ostringstream msg ;
    msg << "Kokkos::DeviceNUMA::" << method << " FAILED: " ;

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
// FULL configuration: use every available CPU
// MOST configuration: use all but one CPU of every NUMA node

void DeviceNUMAInternal::initialize( DeviceNUMA::UtilizationOptions config )
{
  const bool ok_inactive   = 0 == m_thread_count ;
  const bool ok_query_node = 0 != m_cpu_per_node ;

  bool ok_spawn_threads = ok_inactive && ok_query_node ;

  //------------------------------------

  if ( ok_spawn_threads ) {

    switch( config ) {
    case DeviceNUMA::FULL : m_thread_per_node = m_cpu_per_node ; break ;
    case DeviceNUMA::MOST : m_thread_per_node = m_cpu_per_node - 1 ; break ;
    default: break ;
    }

    m_thread_count = m_node_count * m_thread_per_node ;

    //------------------------------------
    // Initialize thread rank and fan-in / fan-out span of threads

    {
      int count = 1 ;

      for ( size_type rank = 0 ; rank < m_thread_count ; ++rank ) {

        m_thread[rank].m_rank      = rank ;
        m_thread[rank].m_fan_begin = m_thread + count ;

        {
          size_type up = 1 ;
          while ( up <= rank )                 { up <<= 1 ; }
          while ( rank + up < m_thread_count ) { up <<= 1 ; ++count ; }
        }

        m_thread[rank].m_fan_end = m_thread + count ;
      }
    }

    //------------------------------------
    // Spawn threads in reverse order for barrier
    {
      m_worker = & m_worker_block ;

      // Create threads last-to-first for start up fan-in barrier

      for ( size_type rank = m_thread_count ; ok_spawn_threads && 1 < rank ; ) {

        DeviceNUMAThread * const thread = m_thread + --rank ;

        thread->set( DeviceNUMAThread::ThreadActive );

        if ( device_numa_thread_spawn( thread ) ) {
          thread->wait( DeviceNUMAThread::ThreadActive );

          ok_spawn_threads =
            DeviceNUMAThread::ThreadInactive == thread->m_state ;
        }
        else {
          thread->set( DeviceNUMAThread::ThreadNull );
          ok_spawn_threads = false ;
        }
      }

      m_worker = NULL ;
    }

    //------------------------------------

    if ( ! ok_spawn_threads ) {
      finalize();
    }
  }

  if ( ! ok_spawn_threads ) {
    std::string msg ;
    msg.append( "Kokkos::DeviceNUMA::initialize() FAILED: " );

    if ( ! ok_inactive ) {
      msg.append( "Device is already active" );
    }
    else if ( ! ok_query_node ) {
      msg.append( "Could not query 'hardware locality (hwloc)'" );
    }
    else {
      msg.append( "Spawning or cpu-binding the threads" );
    }

    throw std::runtime_error( msg );
  }
}

//----------------------------------------------------------------------------

void DeviceNUMAInternal::finalize()
{
  verify_inactive("finalize()");

  while ( 1 < m_thread_count ) {
    DeviceNUMAThread * const thread = m_thread + --m_thread_count ;

    if ( DeviceNUMAThread::ThreadInactive == thread->m_state ) {
      thread->set(  DeviceNUMAThread::ThreadTerminating );
      thread->wait( DeviceNUMAThread::ThreadTerminating );
    }

    thread->m_rank      = 0 ;
    thread->m_fan_begin = NULL ;
    thread->m_fan_end   = NULL ;
  }

  m_thread_count    = 0 ;
  m_thread_per_node = 0 ;
}

//----------------------------------------------------------------------------

void DeviceNUMAInternal::execute( const DeviceNUMAWorker & worker )
{
  verify_inactive("execute(...)");

  m_worker = & worker ;

  DeviceNUMAThread * const thread_beg = m_thread + 1 ;
  DeviceNUMAThread *       thread     = m_thread + m_thread_count ;

  while ( thread_beg < thread ) {
    (--thread)->set( DeviceNUMAThread::ThreadActive );
  }

  // This thread is the root thread of the pool.
  worker.execute_on_thread( m_thread[0] );

  m_worker = NULL ;
}

//----------------------------------------------------------------------------

void DeviceNUMAInternal::block()
{
  verify_inactive("block()");

  device_numa_thread_lock();

  m_worker = & m_worker_block ;

  DeviceNUMAThread * const thread_beg = m_thread + 1 ;
  DeviceNUMAThread *       thread     = m_thread + m_thread_count ;

  while ( thread_beg < thread ) {
    (--thread)->set( DeviceNUMAThread::ThreadActive );
  }

  // All non-root threads are now blocked
}

//----------------------------------------------------------------------------

void DeviceNUMAInternal::unblock()
{
  if ( & m_worker_block != m_worker ) {
    std::ostringstream msg ;
    msg << "Kokkos::DeviceNUMA::unblock() FAILED: " ;
    msg << "Device is not blocked" ;
    throw std::runtime_error( msg.str() );
  }

  device_numa_thread_unlock();

  m_thread->barrier();

  m_worker = NULL ;
}

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

DeviceNUMA::size_type
DeviceNUMAWorker::work_per_thread( DeviceNUMA::size_type work_count )
{
  const DeviceNUMA::size_type
    thread_count = DeviceNUMAInternal::singleton().m_thread_count ;
  return ( work_count + thread_count - 1 ) / thread_count ;
}

//----------------------------------------------------------------------------

} //  namespace Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void DeviceNUMA::block()
{ Impl::DeviceNUMAInternal::singleton().block(); }

void DeviceNUMA::unblock()
{ Impl::DeviceNUMAInternal::singleton().unblock(); }

void DeviceNUMA::execute( const Impl::DeviceNUMAWorker & worker )
{ Impl::DeviceNUMAInternal::singleton().execute( worker ); }

void DeviceNUMA::finalize()
{ Impl::DeviceNUMAInternal::singleton().finalize(); }

void DeviceNUMA::initialize( DeviceNUMA::UtilizationOptions config )
{ Impl::DeviceNUMAInternal::singleton().initialize( config ); }

/*--------------------------------------------------------------------------*/

} // namespace Kokkos


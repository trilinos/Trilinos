/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

/*--------------------------------------------------------------------------*/
/* Kokkos interfaces */

#include <Kokkos_DeviceNUMA.hpp>
#include <impl/Kokkos_MemoryInfo.hpp>
#include <DeviceNUMA/Kokkos_DeviceNUMA_Internal.hpp>

/*--------------------------------------------------------------------------*/
/* Standard 'C' libraries */
#include <stdlib.h>

/* Standard 'C++' libraries */
#include <iostream>
#include <stdexcept>
#include <sstream>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

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

  DeviceNUMAWorkerBlock m_worker_block ;
  size_type             m_node_count ;
  size_type             m_core_per_node ;
  size_type             m_thread_count ;
  size_type             m_thread_per_node ;
  DeviceNUMAThread      m_thread[ THREAD_COUNT_MAX + 1 ];

  const DeviceNUMAWorker * volatile m_worker ;

  ~DeviceNUMAInternal();
  DeviceNUMAInternal();

  void verify_inactive( const char * const method ) const ;

  size_type detect_core_count() const
  { return m_node_count * m_core_per_node ; }

  void initialize( DeviceNUMA::UtilizationStrategy ,
                   DeviceNUMA::size_type manually_set_thread_count );
  void finalize();
  void block();
  void unblock();

  void execute( const DeviceNUMAWorker & worker );

  static DeviceNUMAInternal & singleton();
};

DeviceNUMAInternal & DeviceNUMAInternal::singleton()
{
  static DeviceNUMAInternal self ;
  return self ;
}

DeviceNUMAInternal::~DeviceNUMAInternal()
{}

DeviceNUMAInternal::DeviceNUMAInternal()
  : m_worker_block()
  , m_node_count   ( device_numa_node_count() )
  , m_core_per_node( device_numa_core_per_node() )
  , m_thread_count ( 0 )
  , m_thread_per_node( 0 )
{
  m_worker = NULL ;

  if ( m_node_count && m_core_per_node ) {
    // Bind this master thread to NUMA node 0
    device_numa_bind_this_thread_to_node( 0 );
  }
}

//----------------------------------------------------------------------------
// Driver for each created thread

void DeviceNUMAThread::driver()
{
  try {
    const DeviceNUMAInternal & pool = DeviceNUMAInternal::singleton();

    //------------------------------------
    // this_thread is in the Active state.
    // The master thread is waiting for this_thread to leave the Active state.
    //------------------------------------
    // Migrate myself to the proper NUMA node and run

    if ( device_numa_bind_this_thread_to_node(m_thread_rank/pool.m_thread_per_node) ) {

      while ( ThreadActive == m_state ) {

        // When the work is complete the state will be Inactive or Terminate
        pool.m_worker->execute_on_thread( *this );

        // If this_thread is in the Inactive state then wait for activation.
        wait( ThreadInactive );
      }
    }

    set( ThreadNull );
  }
  catch( const std::exception & x ) {
    std::cerr << "Thread " << m_thread_rank << " uncaught exception : "
              << x.what() << std::endl ;
    std::terminate();
  }
  catch( ... ) {
    std::cerr << "Thread " << m_thread_rank << " uncaught exception"
              << std::endl ;
    std::terminate();
  }
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

void DeviceNUMAInternal::initialize(
  DeviceNUMA::UtilizationStrategy config ,
  DeviceNUMA::size_type manually_set_thread_count )
{
  const bool ok_inactive   = 0 == m_thread_count ;
  const bool ok_query_node = 0 != m_node_count && 0 != m_core_per_node ;

  //------------------------------------

  if ( ok_inactive ) {

    switch( config ) {
    case DeviceNUMA::DETECT_AND_USE_ALL_CORES :
      if ( ok_query_node ) {
        m_thread_per_node = m_core_per_node ;
        m_thread_count    = m_node_count * m_thread_per_node ;
      }
      break ;

    case DeviceNUMA::DETECT_AND_USE_MOST_CORES :
      if ( ok_query_node ) {
        m_thread_per_node = m_core_per_node - 1 ;
        m_thread_count    = m_node_count * m_thread_per_node ;
      }
      break ;

    case DeviceNUMA::MANUALLY_SET_THREAD_COUNT :
      m_thread_count = manually_set_thread_count ;
      if ( ok_query_node ) {
        m_thread_per_node = ( m_thread_count + m_node_count - 1 )
                            / m_node_count ;
      }
      else {
        m_thread_per_node = m_thread_count ;
      }
      break ;
    default: break ;
    }
  }

  bool ok_spawn_threads = ok_inactive && 0 < m_thread_count ;

  //------------------------------------

  if ( ok_spawn_threads ) {

    // Initialize thread rank and fan-in / fan-out span of threads

    {
      int count = 1 ;

      for ( size_type rank = 0 ; rank < m_thread_count ; ++rank ) {

        m_thread[rank].m_thread_count = m_thread_count ;
        m_thread[rank].m_thread_rank  = rank ;
        m_thread[rank].m_thread_reverse_rank = m_thread_count - ( rank + 1 );
        m_thread[rank].m_fan_begin    = m_thread + count ;

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
      msg.append( "Could not query NUMA node count or cores per node" );
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

    thread->m_thread_count        = 0 ;
    thread->m_thread_rank         = 0 ;
    thread->m_thread_reverse_rank = 0 ;
    thread->m_fan_begin           = NULL ;
    thread->m_fan_end             = NULL ;
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
}

void DeviceNUMAInternal::unblock()
{
  if ( & m_worker_block != m_worker ) {
    std::ostringstream msg ;
    msg << "Kokkos::DeviceNUMA::unblock() FAILED: " ;
    msg << "Device is active" ;
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

void DeviceNUMAWorker::execute( const Impl::DeviceNUMAWorker & worker )
{ Impl::DeviceNUMAInternal::singleton().execute( worker ); }

//----------------------------------------------------------------------------

} //  namespace Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void DeviceNUMA::block()
{ Impl::DeviceNUMAInternal::singleton().block(); }

void DeviceNUMA::unblock()
{ Impl::DeviceNUMAInternal::singleton().unblock(); }

void DeviceNUMA::finalize()
{ Impl::DeviceNUMAInternal::singleton().finalize(); }

void DeviceNUMA::initialize( DeviceNUMA::UtilizationStrategy config ,
                             DeviceNUMA::size_type manually_set_thread_count )
{ Impl::DeviceNUMAInternal::singleton().initialize( config , manually_set_thread_count ); }

DeviceNUMA::size_type DeviceNUMA::detect_core_count()
{ return Impl::DeviceNUMAInternal::singleton().detect_core_count(); }

/*--------------------------------------------------------------------------*/

} // namespace Kokkos


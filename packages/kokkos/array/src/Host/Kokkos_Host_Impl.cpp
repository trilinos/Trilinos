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

#include <Kokkos_Host.hpp>
#include <Host/Kokkos_Host_Internal.hpp>

/*--------------------------------------------------------------------------*/
/* Standard 'C' libraries */
#include <stdlib.h>

/* Standard 'C++' libraries */
#include <limits>
#include <iostream>
#include <stdexcept>
#include <sstream>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

HostInternal::~HostInternal()
{}

HostInternal::HostInternal()
  : m_worker_block()
  , m_node_rank( std::numeric_limits<size_type>::max() )
  , m_node_count( 0 )
  , m_node_core_count( 0 )
  , m_page_size( 16 /* default alignment */ )
  , m_thread_count( 1 )
{
  m_worker = NULL ;

  // Master thread:
  m_thread[0].m_thread_count = 1 ;
  m_thread[0].m_thread_rank  = 0 ;
  m_thread[0].m_fan_begin     = NULL ;
  m_thread[0].m_fan_end       = NULL ;
}

//----------------------------------------------------------------------------
// Driver for each created thread

void HostThread::driver()
{
  try {
    const HostInternal & pool = HostInternal::singleton();

    //------------------------------------
    // this_thread is in the Active state.
    // The master thread is waiting for this_thread to leave the Active state.
    //------------------------------------
    // Migrate myself to the proper NUMA node and run

    if ( pool.bind_to_node( *this ) ) {

      while ( ThreadActive == m_state ) {

        // When the work is complete the state will be Inactive or Terminate
        pool.execute( *this );

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

void HostInternal::verify_inactive( const char * const method ) const
{
  if ( NULL != m_worker ) {

    std::ostringstream msg ;
    msg << "Kokkos::Host::" << method << " FAILED: " ;

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
//  Default cannot bind.
//  The 'hwloc' implementation overloads this function

bool HostInternal::bind_to_node( const HostThread & ) const
{ return false ; }

//----------------------------------------------------------------------------

bool HostInternal::spawn_threads()
{
  // If the process is bound to a particular node
  // then only use cores belonging to that node.
  // Otherwise use all nodes and all their cores.

  bool ok_spawn_threads = true ;

  // Initialize thread rank and fan-in / fan-out span of threads

  {
    int count = 1 ;

    for ( size_type rank = 0 ; rank < m_thread_count ; ++rank ) {

      m_thread[rank].m_thread_count = m_thread_count ;
      m_thread[rank].m_thread_rank  = rank ;
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

      HostThread * const thread = m_thread + --rank ;

      thread->set( HostThread::ThreadActive );

      if ( spawn( thread ) ) {
        thread->wait( HostThread::ThreadActive );

        ok_spawn_threads = HostThread::ThreadInactive == thread->m_state ;
      }
      else {
        thread->set( HostThread::ThreadNull );
        ok_spawn_threads = false ;
      }
    }

    m_worker = NULL ;
  }

  //------------------------------------

  if ( ! ok_spawn_threads ) {
    finalize();
  }

  return ok_spawn_threads ;
}

void HostInternal::initialize( const Host::DetectAndUseCores )
{
  if ( m_node_count ) {
    // If pinned to a node then only use cores of that node

    size_t thread_count = m_node_rank < m_node_count
                        ? m_node_core_count : m_node_count * m_node_core_count ;

    initialize( Host::SetThreadCount( thread_count ) );
  }
}

void HostInternal::initialize( const Host::SetThreadCount config )
{
  const bool ok_inactive = 1 == m_thread_count ;
        bool ok_spawn_threads = false ;

  if ( ok_inactive ) {
    m_thread_count = 0 < config.thread_count ? config.thread_count : 1 ;

    ok_spawn_threads = spawn_threads();
  }

  if ( ! ok_spawn_threads ) {
    std::string msg ;
    msg.append( "Kokkos::Host::initialize() FAILED: " );

    if ( ! ok_inactive ) {
      msg.append( "Device is already active" );
    }
    else {
      msg.append( "Spawning or cpu-binding the threads" );
    }

    throw std::runtime_error( msg );
  }

  // Bind this master thread to NUMA node
  bind_to_node( m_thread[0] );
}

//----------------------------------------------------------------------------

void HostInternal::finalize()
{
  verify_inactive("finalize()");

  // Release and clear worker threads:
  while ( 1 < m_thread_count ) {
    HostThread * const thread = m_thread + --m_thread_count ;

    if ( HostThread::ThreadInactive == thread->m_state ) {
      thread->set(  HostThread::ThreadTerminating );
      thread->wait( HostThread::ThreadTerminating );
    }

    thread->m_thread_count = 0 ;
    thread->m_thread_rank  = 0 ;
    thread->m_fan_begin    = NULL ;
    thread->m_fan_end      = NULL ;
  }

  // Reset master thread:
  m_thread->m_thread_count = 1 ;
  m_thread->m_thread_rank  = 0 ;
  m_thread->m_fan_begin    = NULL ;
  m_thread->m_fan_end      = NULL ;

  m_thread_count = 1 ;
}

//----------------------------------------------------------------------------

inline
void HostInternal::execute( const HostThreadWorker<void> & worker )
{
  verify_inactive("execute(...)");

  m_worker = & worker ;

  HostThread::activate( m_thread + 1 , m_thread + m_thread_count );

  // This thread is the root thread of the pool.
  worker.execute_on_thread( m_thread[0] );

  m_worker = NULL ;
}

void HostThreadWorker<void>::execute( const HostThreadWorker<void> & worker )
{ HostInternal::singleton().execute( worker ); }

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

void Host::finalize()
{ Impl::HostInternal::singleton().finalize(); }

void Host::initialize( const Host::DetectAndUseCores config )
{ Impl::HostInternal::singleton().initialize( config ); }

void Host::initialize( const Host::SetThreadCount config )
{ Impl::HostInternal::singleton().initialize( config ); }

Host::size_type Host::detect_core_count()
{
  const Impl::HostInternal & h = Impl::HostInternal::singleton();
  return h.m_node_count * h.m_node_core_count ;
}

Host::size_type Host::detect_node_count()
{ return Impl::HostInternal::singleton().m_node_count ; }

Host::size_type Host::detect_node_core_count()
{ return Impl::HostInternal::singleton().m_node_core_count ; }

Host::size_type Host::detect_memory_page_size()
{ return Impl::HostInternal::singleton().m_page_size ; }

/*--------------------------------------------------------------------------*/

} // namespace Kokkos


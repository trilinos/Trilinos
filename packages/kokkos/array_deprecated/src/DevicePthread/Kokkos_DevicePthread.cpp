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

#include <stdlib.h>
#include <errno.h>
#include <pthread.h>
#include <sched.h>

#include <iostream>
#include <stdexcept>
#include <sstream>

#include <Kokkos_DevicePthread.hpp>
#include <impl/Kokkos_MemoryInfo.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
// Performance critical function: thread waits while control == flag.
// The control value will be changed by a single designated other thread.

void DevicePthreadController::wait( DevicePthreadController::Control flag )
{
  const long value = flag ;
  while ( value == m_control ) {
    sched_yield();
  }
}

//----------------------------------------------------------------------------

class DevicePthreadWorkerBlock : public DevicePthreadWorker {
public:
  void execute_on_thread( DevicePthreadController & ) const ;

  DevicePthreadWorkerBlock() {}
  ~DevicePthreadWorkerBlock() {}

  static DevicePthreadWorkerBlock & singleton();
};

DevicePthreadWorkerBlock & DevicePthreadWorkerBlock::singleton()
{
  static DevicePthreadWorkerBlock self ;
  return self ;
}

//----------------------------------------------------------------------------

class DevicePthreadPool {
public:

  enum { THREAD_COUNT_MAX = 127 };

  typedef DevicePthread::size_type size_type ;

  pthread_mutex_t             m_lock ;
  const DevicePthreadWorker * m_worker ;
  size_type                   m_thread_count ;
  DevicePthreadController     m_thread[ THREAD_COUNT_MAX + 1 ];


  ~DevicePthreadPool() { pthread_mutex_destroy( & m_lock ); }

  DevicePthreadPool()
    {
      pthread_mutex_init( & m_lock , NULL );
      m_worker = NULL ;
      m_thread_count = 0 ;
    }

  void verify_inactive( const char * const method ) const ;

  void initialize( size_type );
  void finalize();
  void block();
  void unblock();

  void execute( const DevicePthreadWorker & worker );

  static void * pthread_driver( void * );

  static DevicePthreadPool & singleton();
};

DevicePthreadPool & DevicePthreadPool::singleton()
{
  static DevicePthreadPool self ;
  return self ;
}

//----------------------------------------------------------------------------
// Driver for each created pthread

void * DevicePthreadPool::pthread_driver( void * arg )
{
  DevicePthreadPool       & pool        = DevicePthreadPool::singleton();
  DevicePthreadController & this_thread = *((DevicePthreadController *) arg);

  do {

    pool.m_worker->execute_on_thread( this_thread );

    // this_thread is now in the DevicePthreadController::Inactive state.
    // Wait until it is reactivated
    this_thread.wait( DevicePthreadController::Inactive );

  } while ( pool.m_worker );

  // If activated with no worker then terminate.

  this_thread.barrier();

  return NULL ;
}

//----------------------------------------------------------------------------

void DevicePthreadPool::verify_inactive( const char * const method ) const
{
  if ( NULL != m_worker ) {
    DevicePthreadWorkerBlock & worker =
      DevicePthreadWorkerBlock::singleton();

    std::ostringstream msg ;
    msg << "Kokkos::DevicePthread::" << method << " FAILED: " ;

    if ( & worker == m_worker ) {
      msg << "Device is blocked" ;
    }
    else {
      msg << "Functor is executing." ;
    }
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------------

void DevicePthreadPool::initialize( size_type nthreads )
{
  const bool error_active = m_thread_count ;
  const bool error_size   = nthreads < 1 || THREAD_COUNT_MAX < nthreads ;

  bool error_pthread = false ;

  if ( ! error_active && ! error_size ) {

    pthread_attr_t attr ;
  
    error_pthread =
      pthread_attr_init( & attr ) ||
      pthread_attr_setscope(       & attr, PTHREAD_SCOPE_SYSTEM ) ||
      pthread_attr_setdetachstate( & attr, PTHREAD_CREATE_DETACHED );

    if ( ! error_pthread ) {
      DevicePthreadWorker noop_worker ;

      size_type thread_rank = 0 ;
      int count = 1 ;

      /* Initialize threads with fan-in / fan-out span of threads */

      for ( thread_rank = 0 ; thread_rank <= nthreads ; ++thread_rank ) {
        DevicePthreadController * const thread = m_thread + thread_rank ;

        thread->m_thread_fan = m_thread + count ;
        thread->m_reduce     = NULL ;
        thread->m_rank       = thread_rank ;
        thread->m_control    = DevicePthreadController::Active ;

        {
          size_type up = 1 ;
          while ( up <= thread_rank )           { up <<= 1 ; }
          while ( thread_rank + up < nthreads ) { up <<= 1 ; ++count ; }
        }
      }

      m_thread_count = nthreads ;
      m_worker       = & noop_worker ;

      /* Create threads last-to-first for start up fan-in barrier */

      for ( thread_rank = nthreads ; ! error_pthread && 1 < thread_rank ; ) {
        DevicePthreadController * const thread = m_thread + --thread_rank ;

        pthread_t pt ;

        if ( pthread_create( & pt, & attr, & pthread_driver, thread ) ) {
          // Spawn failed, set inactive
          thread->set( DevicePthreadController::Inactive );
          error_pthread = true ;
        }
        else {
          thread->wait( DevicePthreadController::Active ); // Wait for self-deactivation
        }
      }

      m_worker = NULL ;

      if ( ! error_pthread ) {
        // All threads successfully spawned
        noop_worker.execute_on_thread( m_thread[0] );
      }
      else {
        /* If any thread-spawn failed, terminate the created threads */
        while ( thread_rank < --( m_thread_count ) ) {
          DevicePthreadController * const thread = m_thread + m_thread_count ;
          thread->set(  DevicePthreadController::Active ); // Re-activate
          thread->wait( DevicePthreadController::Active ); // Wait for termination
        }
        m_thread_count = 0 ;
      }

      pthread_attr_destroy( & attr );
    }
  }

  if ( error_active || error_size || error_pthread ) {
    std::ostringstream msg ;
    msg << "Kokkos::DevicePthread::initialize(" << nthreads << ") FAILED: " ;
    if ( error_active ) {
      msg << "Device is already active" ;
    }
    else if ( error_size ) {
      msg << "Thread count out of range [1.." << THREAD_COUNT_MAX << ")" ;
    }
    else {
      msg << "Error returned from a pthreads system function" ;
    }
  }
}

//----------------------------------------------------------------------------

void DevicePthreadPool::finalize()
{
  verify_inactive("finalize()");

  // Activate threads with NULL worker.

  DevicePthreadController * const thread_beg = m_thread + 1 ;
  DevicePthreadController *       thread     = m_thread + m_thread_count ;

  while ( thread_beg < thread ) {
    (--thread)->set( DevicePthreadController::Active );
  }

  m_thread->barrier();
  m_thread_count = 0 ;
}

//----------------------------------------------------------------------------

void DevicePthreadPool::execute( const DevicePthreadWorker & worker )
{
  verify_inactive("execute(...)");

  m_worker = & worker ;

  DevicePthreadController * const thread_beg = m_thread + 1 ;
  DevicePthreadController *       thread     = m_thread + m_thread_count ;

  while ( thread_beg < thread ) {
    (--thread)->set( DevicePthreadController::Active );
  }

  // This thread is the root thread of the pool.
  worker.execute_on_thread( m_thread[0] );

  m_worker = NULL ;
}

//----------------------------------------------------------------------------

void DevicePthreadPool::block()
{
  DevicePthreadWorkerBlock & worker = DevicePthreadWorkerBlock::singleton();

  verify_inactive("block()");

  m_worker = & worker ;

  pthread_mutex_lock( & m_lock );

  m_worker = & worker ;

  DevicePthreadController * const thread_beg = m_thread + 1 ;
  DevicePthreadController *       thread     = m_thread + m_thread_count ;

  while ( thread_beg < thread ) {
    (--thread)->set( DevicePthreadController::Active );
  }

  // All non-root threads are now blocked on a mutex
}

//----------------------------------------------------------------------------

void DevicePthreadPool::unblock()
{
  DevicePthreadWorkerBlock & worker = DevicePthreadWorkerBlock::singleton();

  if ( & worker != m_worker ) {
    std::ostringstream msg ;
    msg << "Kokkos::DevicePthread::unblock() FAILED: " ;
    msg << "Device is not blocked" ;
    throw std::runtime_error( msg.str() );
  }

  pthread_mutex_unlock( & m_lock );

  m_thread->barrier();

  m_worker = NULL ;
}

//----------------------------------------------------------------------------

void DevicePthreadWorkerBlock::execute_on_thread(
  DevicePthreadController & this_thread ) const
{
  DevicePthreadPool & pool = DevicePthreadPool::singleton();

  pthread_mutex_lock(   & pool.m_lock );
  pthread_mutex_unlock( & pool.m_lock );

  this_thread.barrier();
}

//----------------------------------------------------------------------------

DevicePthread::size_type
DevicePthreadWorker::work_per_thread( DevicePthread::size_type work_count )
{
  const DevicePthread::size_type
    thread_count = DevicePthreadPool::singleton().m_thread_count ;
  return ( work_count + thread_count - 1 ) / thread_count ;
}

// Default driver performs the synchronization barrier
void DevicePthreadWorker::execute_on_thread(
  DevicePthreadController & this_thread ) const
{ this_thread.barrier(); }

//----------------------------------------------------------------------------

} //  namespace Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void DevicePthread::block()
{ Impl::DevicePthreadPool::singleton().block(); }

void DevicePthread::unblock()
{ Impl::DevicePthreadPool::singleton().unblock(); }

void DevicePthread::execute( const Impl::DevicePthreadWorker & worker )
{ Impl::DevicePthreadPool::singleton().execute( worker ); }

void DevicePthread::finalize()
{ Impl::DevicePthreadPool::singleton().finalize(); }

void DevicePthread::initialize( DevicePthread::size_type nthreads )
{ Impl::DevicePthreadPool::singleton().initialize( nthreads ); }

/*--------------------------------------------------------------------------*/

} // namespace Kokkos


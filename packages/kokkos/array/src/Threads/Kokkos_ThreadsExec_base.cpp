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

#include <KokkosArray_config.h>

#include <Kokkos_Threads.hpp>

/*--------------------------------------------------------------------------*/

#if defined( KOKKOS_HAVE_PTHREAD )

/* Standard 'C' Linux libraries */

#include <pthread.h>
#include <sched.h>
#include <errno.h>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
// Driver for each created pthread

namespace {

pthread_mutex_t host_internal_pthread_mutex = PTHREAD_MUTEX_INITIALIZER ;

struct HostInternalPthreadDriver {
  int (* m_driver )( void * );
  void * m_arg ;
  int    m_flag ;
};

void * host_internal_pthread_driver( void * global )
{
  int ( * const driver )( void * ) = ((HostInternalPthreadDriver*) global )->m_driver ;
  void * const arg                 = ((HostInternalPthreadDriver*) global )->m_arg ;

  // Inform launching process that the data has been extracted.
  ((HostInternalPthreadDriver*) global )->m_flag = 1 ;

  global = 0 ; // Don't touch data again

  (*driver)(arg);

  return NULL ;
}

} // namespace

//----------------------------------------------------------------------------

void ThreadsExec::global_lock()
{
  pthread_mutex_lock( & host_internal_pthread_mutex );
}

void ThreadsExec::global_unlock()
{
  pthread_mutex_unlock( & host_internal_pthread_mutex );
}

// Spawn a thread

bool ThreadsExec::spawn( int (*driver)(void*) , void * arg )
{
  bool result = false ;

  pthread_attr_t attr ;

  if ( 0 == pthread_attr_init( & attr ) ||
       0 == pthread_attr_setscope(       & attr, PTHREAD_SCOPE_SYSTEM ) ||
       0 == pthread_attr_setdetachstate( & attr, PTHREAD_CREATE_DETACHED ) ) {

    pthread_t pt ;

    HostInternalPthreadDriver global = { driver , arg , 0 };

    result = 0 == pthread_create( & pt, & attr, host_internal_pthread_driver, & global );

    // Wait until thread has spawned and used the data
    wait_yield( global.m_flag , 0 );
  }

  pthread_attr_destroy( & attr );

  return result ;
}

//----------------------------------------------------------------------------

bool ThreadsExec::is_process()
{
  static const pthread_t master_pid = pthread_self();

  return pthread_equal( master_pid , pthread_self() );
}

//----------------------------------------------------------------------------

namespace {

template< unsigned N > inline void noop_cycle();

template<> inline void noop_cycle<0>() {}
template< unsigned N > inline void noop_cycle()
{
#if !defined ( KOKKOS_DISABLE_ASM ) && \
    ( defined( __GNUC__ ) || \
      defined( __GNUG__ ) || \
      defined( __INTEL_COMPILER__ ) )

  asm volatile("nop");
  noop_cycle<N-1>();

#else
  sched_yield();
#endif
}

}

void ThreadsExec::wait( volatile int & flag , const int value )
{
  // Issue 'NCycle' no-op operations between checks for the flag to change value.
  enum { NCycle = 1 };
  while ( value == flag ) { noop_cycle< NCycle >(); }
}

void ThreadsExec::wait_yield( volatile int & flag , const int value )
{
  while ( value == flag ) { sched_yield(); }
}

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#elif defined( KOKKOS_HAVE_WINTHREAD )

/* Windows libraries */
#include <windows.h>
#include <process.h>

//----------------------------------------------------------------------------
// Driver for each created pthread

namespace KokkosArray {
namespace Impl {
namespace {

unsigned WINAPI host_internal_winthread_driver( void * arg )
{
  HostInternal & pool = HostInternal::singleton();

  pool.driver( reinterpret_cast<size_t>( arg ) );

  return 0 ;
}

class ThreadLockWindows {
private:
  CRITICAL_SECTION  m_handle ;

  ~ThreadLockWindows()
  { DeleteCriticalSection( & m_handle ); }

  ThreadLockWindows();
  { InitializeCriticalSection( & m_handle ); }

  ThreadLockWindows( const ThreadLockWindows & );
  ThreadLockWindows & operator = ( const ThreadLockWindows & );

public:

  static ThreadLockWindows & singleton();

  void lock()
  { EnterCriticalSection( & m_handle ); }

  void unlock()
  { LeaveCriticalSection( & m_handle ); }
};

ThreadLockWindows & ThreadLockWindows::singleton()
{ static ThreadLockWindows self ; return self ; }

} // namespace <>
} // namespace KokkosArray
} // namespace Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

// Spawn this thread

bool HostInternal::spawn( const size_t thread_rank )
{
  void * const arg = reinterpret_cast<void*>( thread_rank );

  unsigned Win32ThreadID = 0 ;

  HANDLE handle =
    _beginthreadex(0,0,host_internal_winthread_driver,arg,0, & Win32ThreadID );

  return ! handle ;
}

//----------------------------------------------------------------------------

void HostWorkerBlock::execute_on_thread( HostThread & this_thread ) const
{
  ThreadLockWindows & lock = ThreadLockWindows::singleton();
  lock.lock();
  lock.unlock();

  this_thread.barrier();
}

//----------------------------------------------------------------------------
// Performance critical function: thread waits while value == *state

void HostMulticoreThread::wait( const HostMulticoreThread::State flag )
{
  const long value = flag ;
  while ( value == m_state ) {
    Sleep(0);
  }
}

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

bool Host::sleep()
{
  Impl::HostInternal & h = Impl::HostInternal::singleton();
  ThreadLockWindows & lock = ThreadLockWindows::singleton();

  const bool is_ready   = NULL == h.m_worker ;
        bool is_blocked = & h.m_worker_block == h.m_worker ;

  if ( is_ready ) {
    ThreadLockWindows::singleton().lock();

    h.m_worker = & h.m_worker_block ;

    Impl::HostThread::activate( h.m_thread + 1 ,
                                h.m_thread + h.m_thread_count );

    is_blocked = true ;
  }

  return is_blocked ;
}

bool Host::wake()
{
  Impl::HostInternal & h = Impl::HostInternal::singleton();

  const bool is_blocked = & h.m_worker_block != h.m_worker ;
        bool is_ready   = NULL == h.m_worker ;

  if ( is_blocked ) {
    ThreadLockWindows::singleton().unlock();

    h.m_thread->barrier();

    h.m_worker = NULL ;

    is_ready = true ;
  }

  return is_ready ;
}

} // namespace KokkosArray

#else /* NO Threads */

namespace KokkosArray {
namespace Impl {

bool host_thread_is_master() { return true ; }
bool host_thread_spawn( int (*)(void*) , void * ) { return false ; }
void host_thread_wait( volatile int * const , const int ) {}
void host_thread_wait_yield( volatile int * const , const int ) {}
void host_thread_lock() {}
void host_thread_unlock() {}

} // namespace Impl
} // namespace KokkosArray

#endif /* End thread model */


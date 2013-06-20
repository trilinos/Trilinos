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

#include <Host/KokkosArray_Host_Thread.hpp>

#include <KokkosArray_Host.hpp>
#include <Host/KokkosArray_Host_Internal.hpp>

#include <limits>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <iostream>

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

void * host_allocate_not_thread_safe(
  const std::string    & label ,
  const std::type_info & scalar_type ,
  const size_t           scalar_size ,
  const size_t           scalar_count );

void host_decrement_not_thread_safe( const void * ptr );

//----------------------------------------------------------------------------

namespace {

struct HostThreadSentinel {
  HostThreadSentinel() {}
  ~HostThreadSentinel();
};

HostThreadSentinel::~HostThreadSentinel()
{
  unsigned i = 0 ;
  for ( ; i < HostThread::max_thread_count && 0 == HostThread::get_thread(i) ; ++i ) {}
  if ( i < HostThread::max_thread_count ) {
    std::cerr << "KokkosArray::Impl::HostThread WARNING : "
              << "existing with live HostThread objects"
              << std::endl ;
  }
}

}

//----------------------------------------------------------------------------

HostThread * HostThread::m_thread[ HostThread::max_thread_count ];
int          HostThread::m_relations = 0 ;

//----------------------------------------------------------------------------

void HostThread::warn_destroy_with_reduce()
{
  std::cerr << "KokkosArray::Impl::HostThread WARNING : destroyed with allocated reduction memory"
            << std::endl ;
}

void HostThread::set_topology( const unsigned thread_rank , const unsigned thread_count ,
                               const unsigned gang_rank ,   const unsigned gang_count ,
                               const unsigned worker_rank , const unsigned worker_count )
{
  if ( m_thread[ thread_rank ] != this ) {
    std::ostringstream msg ;
    msg << "KokkosArray::Impl::HostThread::set_topology("
        << thread_rank << "," << thread_count << ","
        << gang_rank << "," << gang_count << ","
        << worker_rank << "," << worker_count << ")"
        << " ERROR : Does not match previously set thread rank" ;
    if ( m_thread[ m_thread_rank ] ) {
     msg << " : ("
         << m_thread[ thread_rank ]->m_thread_rank << ","
         << m_thread[ thread_rank ]->m_thread_count << ","
         << m_thread[ thread_rank ]->m_gang_rank << ","
         << m_thread[ thread_rank ]->m_gang_count << ","
         << m_thread[ thread_rank ]->m_worker_rank << ","
         << m_thread[ thread_rank ]->m_worker_count << ")" ;
    }
    throw std::runtime_error( msg.str() );
  }

  m_thread_rank  = thread_rank ;
  m_thread_count = thread_count ;

  m_gang_rank    = gang_rank ;
  m_gang_count   = gang_count ;

  m_worker_rank  = worker_rank ;
  m_worker_count = worker_count ;
}


void HostThread::resize_reduce( unsigned size )
{
  if ( m_reduce ) {
    host_decrement_not_thread_safe( m_reduce );
    m_reduce = 0 ;
  }

  if ( size ) {

    m_reduce = host_allocate_not_thread_safe( "reduce_scratch_space" , typeid(unsigned char) , 1 , size );

    // Guaranteed multiple of 'unsigned'

    unsigned * ptr = (unsigned *)( m_reduce );
    unsigned * const end = ptr + size / sizeof(unsigned);

    // touch on this thread
    while ( ptr < end ) *ptr++ = 0 ;
  }
}

void HostThread::set_thread_relationships()
{
  if ( m_relations ) return ;

  bool error_gaps   = false ;
  bool error_thread = false ;

  const unsigned gang_count   = HostThread::m_thread[0] ? HostThread::m_thread[0]->m_gang_count : 0 ;
  const unsigned thread_count = HostThread::m_thread[0] ? HostThread::m_thread[0]->m_thread_count : 0 ;

  {
    unsigned i = 0 ;
    unsigned j = 0 ;

    // Non-zero entries:
    for ( i = 0 ; i < HostThread::max_thread_count && 0 != HostThread::m_thread[i]; ++i ) {}

    // Followed by zero entries:
    for ( j = i ; j < HostThread::max_thread_count && 0 == HostThread::m_thread[j] ; ++j ) {}

    if ( 0 == i && HostThread::max_thread_count == j ) return ;

    error_gaps   = j < HostThread::max_thread_count ;
    error_thread = i != thread_count ;
  }

  if ( ! error_gaps && ! error_thread ) {

    for ( unsigned g = 0 , r = 0 ; g < gang_count && r < thread_count ; ++g ) {

      const unsigned worker_count = HostThread::m_thread[r]->m_worker_count ;

      for ( unsigned w = 0 ; w < worker_count ; ++w , ++r ) {
        if ( r != HostThread::m_thread[r]->m_thread_rank ||
             g != HostThread::m_thread[r]->m_gang_rank ||
             w != HostThread::m_thread[r]->m_worker_rank ||
             thread_count != HostThread::m_thread[r]->m_thread_count ||
             gang_count   != HostThread::m_thread[r]->m_gang_count ||
             worker_count != HostThread::m_thread[r]->m_worker_count ) {
          error_thread = true ;
        }
      }
    }
  }
  else {
    std::ostringstream msg ;
    msg << "KokkosArray::Impl::HostThread::set_thread_ralationships ERROR :" ;
    if ( error_gaps ) {
      msg << " Gaps in thread array" ;
    }
    else {
      msg << " Inconsistent thread topology" ;
    }
    throw std::runtime_error( msg.str() );
  }

  // Assign fan in/out for each thread
  for ( unsigned i = 0 ; i < thread_count ; ++i ) {

    HostThread & thread = * HostThread::m_thread[i] ; // i == thread.m_thread_rank

    unsigned fan_count = 0 ;

    // Intra-gang reduction:
    for ( unsigned n = 1 ; thread.m_worker_rank + n < thread.m_worker_count ; n <<= 1 ) {

      if ( n & thread.m_worker_rank ) break ;

      thread.m_fan[ fan_count++ ] = HostThread::m_thread[ thread.m_thread_rank + n ];
    }

    // Inter-gang reduction:
    if ( thread.m_worker_rank == 0 ) {

      for ( unsigned n = 1 ; thread.m_gang_rank + n < thread.m_gang_count ; n <<= 1 ) {

        if ( n & thread.m_gang_rank ) break ;

        // Find thread with:
        //   gank_rank   == thread.m_gang_rank + n
        //   worker_rank == 0
        HostThread ** th = HostThread::m_thread + i ;
        while ( (*th)->m_gang_rank < thread.m_gang_rank + n ) {
          th += (*th)->m_worker_count ;
        }

        thread.m_fan[ fan_count++ ] = *th ;
      }
    }

    thread.m_fan_count = fan_count ;
  }

  m_relations = 1 ;

#if 0

  for ( unsigned i = 0 ; i < HostThread::max_thread_count ; ++i ) {
    if ( 0 != HostThread::m_thread[i] ) {
      HostThread & thread = * HostThread::m_thread[i] ;

      std::cout << "HostThread[" << i << "] :"
                << " rank[ " << thread.m_thread_rank
                << " / " << thread.m_thread_count << " ]"
                << " gang[ " << thread.m_gang_tag
                << " : " << thread.m_gang_rank
                << " / " << thread.m_gang_count << " ]"
                << " worker[ " << thread.m_worker_tag
                << " : " << thread.m_worker_rank
                << " / " << thread.m_worker_count << " ]"
                ;

      if ( thread.m_fan_count ) {
        std::cout << " fan_ranks[" ;
        for ( unsigned j = 0 ; j < thread.m_fan_count ; ++j ) {
          std::cout << " " << thread.m_fan[j]->m_thread_rank ;
        }
        std::cout << " ]" ;
      }

      std::cout << std::endl ;
    }
  }

#endif

}

void HostThread::clear_thread_relationships()
{
  if ( m_relations ) {
    for ( unsigned i = 0 ; i < HostThread::max_thread_count ; ++i ) {
      if ( HostThread::m_thread[i] ) {

        HostThread & thread = * HostThread::m_thread[i] ;

        // Reset to view self as the only thread:

        thread.m_fan_count    = 0 ;
        thread.m_thread_rank  = 0 ;
        thread.m_thread_count = 1 ;
        thread.m_gang_rank    = 0 ;
        thread.m_gang_count   = 1 ;
        thread.m_worker_rank  = 0 ;
        thread.m_worker_count = 1 ;

        for ( unsigned j = 0 ; j < HostThread::max_fan_count ; ++j ) {
          thread.m_fan[j] = 0 ;
        }
      }
    }
    m_relations = 0 ;
  }
}

void HostThread::set_thread( const unsigned rank , HostThread * t )
{
  static const HostThreadSentinel sentinel ;

  HostThread::clear_thread_relationships();

  const bool ok_rank = rank < HostThread::max_thread_count ;
  const bool ok_zero = ok_rank && ( 0 == HostThread::m_thread[ rank ] );

  if ( ok_rank && ok_zero ) {
    HostThread::m_thread[ rank ] = t ;
  }
  else {
    std::ostringstream msg ;
    msg << "KokkosArray::Impl::HostThread::set_thread( "
        << rank << " , ... ) ERROR: " ;
    if ( ! ok_rank ) { msg << " OUT OF BOUNDS" ; }
    else if ( ! ok_zero ) { msg << " ALREADY SET" ; }
    throw std::runtime_error( msg.str() );
  }
}

HostThread * HostThread::clear_thread( const unsigned rank )
{
  HostThread::clear_thread_relationships();

  HostThread * th = 0 ;

  if ( rank < HostThread::max_thread_count ) {
    th = HostThread::m_thread[ rank ] ;
    HostThread::m_thread[ rank ] = 0 ;
  }
  else {
    std::ostringstream msg ;
    msg << "KokkosArray::Impl::HostThread::clear_thread( "
        << rank << " , ... ) ERROR:  OUT OF BOUNDS" ;
    throw std::runtime_error( msg.str() );
  }

  return th ;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#if defined( KOKKOSARRAY_HAVE_PTHREAD )

/*--------------------------------------------------------------------------*/
/* Standard 'C' Linux libraries */
#include <pthread.h>
#include <sched.h>
#include <errno.h>

/*--------------------------------------------------------------------------*/

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------
// Driver for each created pthread

namespace {

pthread_mutex_t host_internal_pthread_mutex = PTHREAD_MUTEX_INITIALIZER ;

void * host_internal_pthread_driver( void * )
{
  try {
    host_thread_driver();
  }
  catch( const std::exception & x ) {
    // mfh 29 May 2012: Doesn't calling std::terminate() seem a
    // little violent?  On the other hand, C++ doesn't define how
    // to transport exceptions between threads (until C++11).
    // Since this is a worker thread, it would be hard to tell the
    // master thread what happened.
    std::cerr << "Worker thread uncaught exception : " << x.what() << std::endl ;
    std::terminate();
  }
  catch( ... ) {
    // mfh 29 May 2012: See note above on std::terminate().
    std::cerr << "Worker thread uncaught exception" << std::endl ;
    std::terminate();
  }

  return NULL ;
}

}

//----------------------------------------------------------------------------

bool host_thread_is_master()
{
  static const pthread_t master_pid = pthread_self();

  return pthread_equal( master_pid , pthread_self() );
}

//----------------------------------------------------------------------------
// Spawn this thread

bool host_thread_spawn()
{
  bool result = false ;

  pthread_attr_t attr ;

  if ( 0 == pthread_attr_init( & attr ) ||
       0 == pthread_attr_setscope(       & attr, PTHREAD_SCOPE_SYSTEM ) ||
       0 == pthread_attr_setdetachstate( & attr, PTHREAD_CREATE_DETACHED ) ) {

    pthread_t pt ;

    result =
      0 == pthread_create( & pt, & attr, host_internal_pthread_driver, 0 );
  }

  pthread_attr_destroy( & attr );

  return result ;
}

//----------------------------------------------------------------------------

namespace {

template< unsigned N > inline void noop_cycle();

template<> inline void noop_cycle<0>() {}
template< unsigned N > inline void noop_cycle()
{
#if defined( __GNUC__ ) || \
    defined( __GNUG__ ) || \
    defined( __INTEL_COMPILER__ )

  asm volatile("nop");
  noop_cycle<N-1>();

#else
  sched_yield();
#endif
}

}

void host_thread_wait( volatile int * const flag , const int value )
{
  // Issue 'NCycle' no-op operations between checks for the flag to change value.
  enum { NCycle = 1 };
  while ( value == *flag ) { noop_cycle< NCycle >(); }
}

void host_thread_wait_yield( volatile int * const flag , const int value )
{
  while ( value == *flag ) { sched_yield(); }
}

void host_thread_lock()
{
  pthread_mutex_lock( & host_internal_pthread_mutex );
}

void host_thread_unlock()
{
  pthread_mutex_unlock( & host_internal_pthread_mutex );
}

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#elif defined( KOKKOSARRAY_HAVE_WINTHREAD )

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
bool host_thread_spawn()     { return false ; }
void host_thread_wait( volatile int * const , const int ) {}
void host_thread_wait_yield( volatile int * const , const int ) {}
void host_thread_lock() {}
void host_thread_unlock() {}

} // namespace Impl
} // namespace KokkosArray

#endif /* End thread model */


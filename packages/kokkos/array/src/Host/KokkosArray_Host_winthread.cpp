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

/*--------------------------------------------------------------------------*/
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


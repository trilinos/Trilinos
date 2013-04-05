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

void * host_internal_pthread_driver( void * arg )
{
  HostInternal::singleton().driver( reinterpret_cast<size_t>( arg ) );

  return NULL ;
}

pthread_mutex_t host_internal_pthread_mutex = PTHREAD_MUTEX_INITIALIZER ;

}

//----------------------------------------------------------------------------

bool HostInternal::is_master_thread() const
{
  static const pthread_t master_pid = pthread_self();

  return pthread_equal( master_pid , pthread_self() );
}

//----------------------------------------------------------------------------
// Spawn this thread

bool HostInternal::spawn( const size_t thread_rank )
{
  bool result = false ;

  pthread_attr_t attr ;
  
  if ( 0 == pthread_attr_init( & attr ) ||
       0 == pthread_attr_setscope(       & attr, PTHREAD_SCOPE_SYSTEM ) ||
       0 == pthread_attr_setdetachstate( & attr, PTHREAD_CREATE_DETACHED ) ) {

    void * const arg = reinterpret_cast<void*>( thread_rank );

    pthread_t pt ;

    result =
      0 == pthread_create( & pt, & attr, host_internal_pthread_driver, arg );
  }

  pthread_attr_destroy( & attr );

  return result ;
}

//----------------------------------------------------------------------------

void host_thread_wait( volatile int * const flag , const int value )
{
  while ( value == *flag ) sched_yield();
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


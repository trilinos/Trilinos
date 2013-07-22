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

#include <iostream>
#include <Kokkos_Threads.hpp>
#include <KokkosArray_hwloc.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

void * host_allocate_not_thread_safe(
  const std::string    & label ,
  const std::type_info & scalar_type ,
  const size_t           scalar_size ,
  const size_t           scalar_count );

void host_decrement_not_thread_safe( const void * ptr );

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
namespace {

enum { CURRENT_MEMORY_FUNCTOR_SIZE = 0x08000 /* 32k words , 128k bytes */ };

int s_current_functor_memory[ CURRENT_MEMORY_FUNCTOR_SIZE ];

ThreadsExec                  s_threads_process ;
ThreadsExec                * s_threads_exec[  ThreadsExec::MAX_THREAD_COUNT ];
std::pair<unsigned,unsigned> s_threads_coord[ ThreadsExec::MAX_THREAD_COUNT ];
std::string                  s_exception_msg ;

unsigned s_threads_count       = 0 ;
unsigned s_threads_reduce_size = 0 ;
unsigned s_threads_shared_size = 0 ;

void (* volatile s_current_function)( Threads , const void * );

struct Sentinel {
  Sentinel() {}
  ~Sentinel()
  {
    if ( s_threads_count ||
         s_threads_reduce_size ||
         s_threads_shared_size ||
         s_current_function ||
         s_threads_exec[0] ) {
      std::cerr << "ERROR : Process exiting without calling Kokkos::Threads::terminate()" << std::endl ;
    }
  }
};

} // namespace
} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// C11 threads compatible function:

int ThreadsExec::driver( void * arg )
{
  // If hardware locality library unavailable then pass in the rank.

  const unsigned thread_rank = 
    arg ? *((unsigned*)arg)
        : KokkosArray::hwloc::bind_this_thread( s_threads_count , s_threads_coord );

  if ( s_threads_count <= thread_rank || 0 != s_threads_exec[ thread_rank ] ) {

    // An error occured. Inform process that thread is terminating
    s_threads_process.m_state = ThreadsExec::Terminating ;

    return 0 ;
  }

  {
    ThreadsExec this_thread ;

    this_thread.m_state = ThreadsExec::Active ;

    s_threads_exec[ thread_rank ] = & this_thread ;

    // Inform spawning process that the threads_exec entry has been set.
    s_threads_process.m_state = ThreadsExec::Active ;

    while ( ThreadsExec::Active == this_thread.m_state ) {

      try {
        // Call work function
        (*s_current_function)( Threads( this_thread ) , s_current_functor_memory );
      }
      catch( const std::exception & x ) {
        std::ostringstream msg ;
        msg << "Kokkos::Threads[" << thread_rank << "] Uncaught exeception : " << x.what() << std::endl ;
        s_exception_msg.append( msg.str() );
      }
      catch( ... ) {
        std::ostringstream msg ;
        msg << "Kokkos::Threads[" << thread_rank << "] Uncaught exeception"  << std::endl ;
        s_exception_msg.append( msg.str() );
      }

      // Deactivate thread and wait for reactivation
      this_thread.m_state = ThreadsExec::Inactive ;
      wait_yield( this_thread.m_state , ThreadsExec::Inactive );
    }

    s_threads_process.m_state = ThreadsExec::Terminating ;

    s_threads_exec[ thread_rank ] = 0 ;
  }

  return 0 ;
}

void execute_function_noop( Threads , const void * ) {}

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

ThreadsExec::~ThreadsExec()
{
  m_reduce      = 0 ;
  m_shared      = 0 ;
  m_state       = ThreadsExec::Terminating ;
  m_fan_size    = 0 ;
  m_team_rank   = 0 ;
  m_team_size   = 0 ;
  m_league_rank = 0 ;
  m_league_rank = 0 ;
  m_thread_size = 0 ;
  m_thread_size = 0 ;

  for ( unsigned i = 0 ; i < MAX_FAN_COUNT ; ++i ) { m_fan[i] = 0 ; }
}

ThreadsExec::ThreadsExec()
  : m_reduce(0)
  , m_shared(0)
  , m_state( ThreadsExec::Terminating )
  , m_fan_size(0)
  , m_team_rank(0)
  , m_team_size(0)
  , m_league_rank(0)
  , m_league_size(0)
  , m_thread_rank(0)
  , m_thread_size(0)
{
  for ( unsigned i = 0 ; i < MAX_FAN_COUNT ; ++i ) { m_fan[i] = 0 ; }
}

void ThreadsExec::set_threads_relationships(
  const std::pair<unsigned,unsigned> team_topo ,
  ThreadsExec * threads[] )
{
  const unsigned thread_count = team_topo.first * team_topo.second ;

  for ( unsigned r = 0 ; r < thread_count ; ++r ) {
    if ( threads[r] == 0 ) {
      KokkosArray::Impl::throw_runtime_exception( std::string("ThreadsExec::set_threads_relationships FAILED : NULL entry" ) );
    }
  }

  for ( unsigned league_rank = 0 , r = 0 ; league_rank < team_topo.first ;  ++league_rank ) {
  for ( unsigned team_rank = 0 ;           team_rank   < team_topo.second ; ++team_rank , ++r ) {

    ThreadsExec & th = * threads[r] ;

    th.m_league_rank = league_rank ;
    th.m_league_size = team_topo.first ;
    th.m_team_rank   = team_rank ;
    th.m_team_size   = team_topo.second ;
    th.m_thread_rank = th.m_team_rank + th.m_team_size * th.m_league_rank ;
    th.m_thread_size = th.m_team_size * th.m_league_size ;

    th.m_fan_size = 0 ;

    // Intra-team reduction:
    for ( int n = 1 ; ( th.m_team_rank + n < th.m_team_size ) &&
                      ( 0 == ( n & th.m_team_rank ) ) ; n <<= 1 ) {
      th.m_fan[ th.m_fan_size++ ] = threads[ ( th.m_team_rank + n ) + ( th.m_league_rank * th.m_team_size ) ];
    }

    // Inter-team (intra-league) reduction:

    if ( th.m_team_rank == 0 ) {

      for ( int n = 1 ; ( th.m_league_rank + n < th.m_league_size ) &&
                        ( 0 == ( n & th.m_league_rank ) ) ; n <<= 1 ) {

        th.m_fan[ th.m_fan_size++ ] = threads[ ( th.m_league_rank + n ) * th.m_team_size ];
      }
    }
  }}
}

void ThreadsExec::execute_sleep( Threads dev , const void * )
{
  ThreadsExec::global_lock();
  ThreadsExec::global_unlock();

  const int n = dev.m_exec.m_fan_size ;

  for ( int i = 0 ; i < n ; ++i ) {
    wait( dev.m_exec.m_fan[i]->m_state , ThreadsExec::Active );
  }

  dev.m_exec.m_state = ThreadsExec::Inactive ;
}

void ThreadsExec::execute_reduce_resize( Threads dev , const void * )
{
  if ( dev.m_exec.m_reduce ) {
    KokkosArray::Impl::host_decrement_not_thread_safe( dev.m_exec.m_reduce );
    dev.m_exec.m_reduce = 0 ;
  }

  if ( s_threads_reduce_size ) {

    dev.m_exec.m_reduce =
      KokkosArray::Impl::host_allocate_not_thread_safe( "reduce_scratch_space" , typeid(unsigned char) , 1 , s_threads_reduce_size );

    // Guaranteed multiple of 'unsigned'

    unsigned * ptr = (unsigned *)( dev.m_exec.m_reduce );
    unsigned * const end = ptr + s_threads_reduce_size / sizeof(unsigned);

    // touch on this thread
    while ( ptr < end ) *ptr++ = 0 ;
  }
}

void ThreadsExec::execute_shared_resize( Threads dev , const void * )
{
  if ( dev.m_exec.m_team_rank ) {
    dev.m_exec.m_shared = 0 ;
  }
  else {

    if ( dev.m_exec.m_shared ) {
      KokkosArray::Impl::host_decrement_not_thread_safe( dev.m_exec.m_shared );
      dev.m_exec.m_shared = 0 ;
    }

    if ( s_threads_shared_size ) {

      dev.m_exec.m_shared =
        KokkosArray::Impl::host_allocate_not_thread_safe( "shared_scratch_space" , typeid(unsigned char) , 1 , s_threads_shared_size );

      // Guaranteed multiple of 'unsigned'

      unsigned * ptr = (unsigned *)( dev.m_exec.m_shared );
      unsigned * const end = ptr + s_threads_shared_size / sizeof(unsigned);

      // touch on this thread
      while ( ptr < end ) *ptr++ = 0 ;
    }
  }
}

}
}

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

void ThreadsExec::verify_is_process( const std::string & name )
{
  if ( ! is_process() ) {
    std::string msg( name );
    msg.append( " FAILED : Called by a worker thread, can only be called by the master process." );
    KokkosArray::Impl::throw_runtime_exception( msg );
  }
}

// Wait for root thread to become inactive
void ThreadsExec::fence()
{
  if ( s_threads_count ) {
    wait_yield( s_threads_exec[0]->m_state , ThreadsExec::Active );

    s_current_function = 0 ;

    if ( s_exception_msg.size() ) {
      KokkosArray::Impl::throw_runtime_exception( s_exception_msg );
    }
  }
}

inline
void ThreadsExec::activate_threads()
{
  for ( unsigned i = s_threads_count ; 0 < --i ; ) {
    s_threads_exec[i]->m_state = ThreadsExec::Active ;
  }
}

/** \brief  Acquire memory for the asynchronous functor */
void * ThreadsExec::acquire( size_t n )
{
  verify_is_process("ThreadsExec::acquire");

  fence();

  if ( CURRENT_MEMORY_FUNCTOR_SIZE < n ) {
    std::ostringstream msg ;
    msg << "ThreadsExec::aquire( " << n << " ) FAILED : requested functor memory larger than "
        << unsigned(CURRENT_MEMORY_FUNCTOR_SIZE) ;
    KokkosArray::Impl::throw_runtime_exception( msg.str() );
  }

  return s_current_functor_memory ;
}

/** \brief  Begin execution of the asynchronous functor */
void ThreadsExec::start( void (*func)( Kokkos::Threads , const void * ) )
{
  verify_is_process("ThreadsExec::start");

  s_exception_msg.clear();

  s_current_function = func ;

  activate_threads();

  if ( s_threads_process.m_thread_size ) {
    (*func)( Threads( s_threads_process ) , s_current_functor_memory );
    s_threads_process.m_state = ThreadsExec::Inactive ;
  }
}

bool ThreadsExec::sleep()
{
  verify_is_process("ThreadsExec::sleep");

  if ( & execute_sleep == s_current_function ) return false ;

  fence();

  ThreadsExec::global_lock();

  s_exception_msg.clear();

  s_current_function = & execute_sleep ;

  activate_threads();
  
  return true ;
}

bool ThreadsExec::wake()
{
  verify_is_process("ThreadsExec::wake");

  if ( & execute_sleep != s_current_function ) return false ;

  ThreadsExec::global_unlock();

  if ( s_threads_process.m_thread_size ) {
    execute_sleep( Threads( s_threads_process ) , s_current_functor_memory );
    s_threads_process.m_state = ThreadsExec::Inactive ;
  }

  fence();

  return true ;
}

//----------------------------------------------------------------------------

void ThreadsExec::execute_serial( void (*func)( Kokkos::Threads , const void * ) )
{
  s_exception_msg.clear();

  s_current_function = func ;

  const unsigned begin = s_threads_process.m_thread_size ? 1 : 0 ;

  for ( unsigned i = s_threads_count ; begin < i ; ) {
    ThreadsExec & th = * s_threads_exec[ --i ];

    th.m_state = ThreadsExec::Active ;

    wait_yield( th.m_state , ThreadsExec::Active );
  }

  if ( s_threads_process.m_thread_size ) {
    s_threads_process.m_state = ThreadsExec::Active ;
    (*func)( Threads( s_threads_process ) , s_current_functor_memory );
    s_threads_process.m_state = ThreadsExec::Inactive ;
  }

  s_current_function = 0 ;
}

//----------------------------------------------------------------------------

void ThreadsExec::resize_reduce_scratch( size_t size )
{
  fence();

  const size_t rem = size % KokkosArray::Impl::MEMORY_ALIGNMENT ;

  if ( rem ) size += KokkosArray::Impl::MEMORY_ALIGNMENT - rem ;

  if ( s_threads_reduce_size < size || 0 == size ) {

    s_threads_reduce_size = size ;

    execute_serial( & execute_reduce_resize );
  }
}

void ThreadsExec::resize_shared_scratch( size_t size )
{
  fence();

  const size_t rem = size % KokkosArray::Impl::MEMORY_ALIGNMENT ;

  if ( rem ) size += KokkosArray::Impl::MEMORY_ALIGNMENT - rem ;

  if ( s_threads_shared_size < size || 0 == size ) {
    s_threads_shared_size = size ;

    execute_serial( & execute_shared_resize );
  }

  for ( unsigned i = 0 ; i < s_threads_count ; ) {
    ThreadsExec & team_th = * s_threads_exec[i] ;

    for ( int j = 0 ; j < team_th.m_team_size ; ++j , ++i ) {
      s_threads_exec[i]->m_shared = team_th.m_shared ;
    }
  }
}

//----------------------------------------------------------------------------

void ThreadsExec::print_configuration( std::ostream & s , const bool detail )
{
  verify_is_process("ThreadsExec::print_configuration");

  fence();

  const std::pair<unsigned,unsigned> core_topo = KokkosArray::hwloc::get_core_topology();
  const unsigned core_size = KokkosArray::hwloc::get_core_capacity();

  s << "Kokkos::Threads hwloc[" << core_topo.first << "x" << core_topo.second << "x" << core_size << "]" ;

  if ( s_threads_exec[0] ) {
    s << " team_league[" << s_threads_exec[0]->m_league_size << "x" << s_threads_exec[0]->m_team_size << "]" ;
    if ( 0 == s_threads_process.m_thread_size ) { s << " Asynchronous" ; }

    if ( detail ) {
      s << " ReduceScratch[" << s_threads_reduce_size << "]"
        << " SharedScratch[" << s_threads_shared_size << "]" ;
      s << std::endl ;
      for ( unsigned i = 0 ; i < s_threads_count ; ++i ) {
        ThreadsExec * const th = s_threads_exec[i] ;
        s << "  Thread" ;
        if ( th ) {
          s << "[" << th->m_league_rank << "." << th->m_team_rank << "]" ;
          if ( th->m_fan_size ) {
            s << " Fan" ;
            for ( int j = 0 ; j < th->m_fan_size ; ++j ) {
              s << " [" << th->m_fan[j]->m_league_rank << "." << th->m_fan[j]->m_team_rank << "]" ;
            }
          }
        }
      }
    }

    s << std::endl ;
  }
  else {
    s << " not initialized" << std::endl ;
  }
}

//----------------------------------------------------------------------------

void ThreadsExec::initialize( 
  const std::pair<unsigned,unsigned> team_topo ,
        std::pair<unsigned,unsigned> core_use )
{
  static const Sentinel sentinel ;

  verify_is_process("ThreadsExec::initialize");

  std::ostringstream msg ;

  msg << "Kokkos::Threads::initialize("
      << " team_topo(" << team_topo.first << "," << team_topo.second << ")"
      << ", core_use(" << core_use.first << "," << core_use.second << ")"
      << " )" ;

  if ( s_threads_count ) {
    msg << " FAILED : Already initialized" ;
    KokkosArray::Impl::throw_runtime_exception( msg.str() );
  }

  //------------------------------------
  // Query hardware topology and capacity, if available.

  const bool                         hwloc_avail  = KokkosArray::hwloc::available();
  const std::pair<unsigned,unsigned> master_coord = KokkosArray::hwloc::get_this_thread_coordinate();
  const std::pair<unsigned,unsigned> core_topo    = KokkosArray::hwloc::get_core_topology();
  const unsigned                     core_cap     = KokkosArray::hwloc::get_core_capacity();
  const unsigned                     capacity     = hwloc_avail ? core_topo.first * core_topo.second * core_cap : 0 ;
  const unsigned                     thread_count = team_topo.first * team_topo.second ;

  //------------------------------------
  // Use HWLOC to determine coordinates for pinning threads.

  if ( hwloc_avail ) {

    if ( capacity         < thread_count ||
         core_topo.first  < core_use.first ||
         core_topo.second < core_use.second ) {
      msg << " FAILED : Requested more cores or threads than HWLOC reports are available "
          << " core_topology(" << core_topo.first << "," << core_topo.second << ")"
          << " thread_capacity(" << capacity << ")" ;
      KokkosArray::Impl::throw_runtime_exception( msg.str() );
    }

    if ( 0 == core_use.first || 0 == core_use.second ) {
      // User requested that we determine best use of cores.

      // Start by assuming use of all available cores
      core_use = core_topo ;

      // Can spawn all requested threads without using a core?
      if ( thread_count <= core_topo.first * ( core_topo.second - 1 ) * core_cap ) {
        --core_use.second ;
      }
    }

    if ( core_use.first < core_topo.first ) {
      // Can omit the use of group of cores and execute work asynchronously

      KokkosArray::Impl::host_thread_mapping( team_topo , core_use , core_topo , s_threads_coord );

      // Don't use master thread's first core coordinate.
      for ( unsigned i = 0 ; i < thread_count ; ++i ) {
        s_threads_coord[i].first = ( s_threads_coord[i].first + master_coord.first + 1 ) % core_topo.first ;
      }
    }
    else if ( core_use.second < core_topo.second ) {
      // Can omit the use of a core and execute work asynchronously

      KokkosArray::Impl::host_thread_mapping( team_topo , core_use , core_topo , s_threads_coord );

      // Don't use the master thread's second core coordinate.
      for ( unsigned i = 0 ; i < thread_count ; ++i ) {
        if ( s_threads_coord[i].first == master_coord.first ) {
          s_threads_coord[i].second = ( s_threads_coord[i].second + master_coord.second + 1 ) % core_topo.second ;
        }
      }
    }
    else {
      // Spawn threads with root thread on the master process' core

      KokkosArray::Impl::host_thread_mapping( team_topo , core_use , core_topo , master_coord , s_threads_coord );
    }

    if ( thread_count == capacity ) {
      // Fully subscribed so claim coordinate #0 for master thread

      s_threads_coord[0] = std::pair<unsigned,unsigned>( ~0u , ~0u );
    }
  }

  //------------------------------------
  // Spawn threads

  {
    const unsigned thread_spawn_begin  = thread_count < capacity ? 0 : 1 ;
    unsigned       thread_spawn_failed = 0 ;

    s_threads_count    = thread_count ;
    s_current_function = & execute_function_noop ; // Initialization work function

    // If not fully utilizing the capacity then spawn threads for asynchronous execution.

    for ( unsigned i = thread_spawn_begin ; i < thread_count ; ++i ) {

      s_threads_process.m_state = ThreadsExec::Inactive ;

      // Spawn thread executing the 'driver' function.
      // If hwloc then will choose its own rank, otherwise specify the rank.
      if ( ThreadsExec::spawn( & ThreadsExec::driver , ( hwloc_avail ? (unsigned*)0 : & i ) ) ) {

        // Thread spawned.  Thread will inform me of its condition as:
        //   Active     = success and has entered its data in the table
        //   Terminate  = failure

        wait_yield( s_threads_process.m_state , ThreadsExec::Inactive );
      }
    }

    // All successfully spawned threads are in the list as [1,thread_count)
    // Wait for them to deactivate before zeroing the worker

    for ( unsigned i = thread_spawn_begin ; i < thread_count ; ++i ) {
      ThreadsExec * const th = s_threads_exec[i] ;
      if ( th ) {
        wait_yield( th->m_state , ThreadsExec::Active );
      }
      else {
        ++thread_spawn_failed ;
      }
    }

    s_current_function = NULL ;

    if ( thread_spawn_failed ) {

      s_threads_count = 0 ;

      msg << " FAILED : with " << thread_spawn_failed << " attempts to spawn threads" ;

      KokkosArray::Impl::throw_runtime_exception( msg.str() );
    }

    if ( thread_spawn_begin ) { // Include the master thread
      KokkosArray::hwloc::bind_this_thread( master_coord );
      s_threads_process.m_state = ThreadsExec::Inactive ;
      s_threads_exec[0] = & s_threads_process ;
    }
  }

  //------------------------------------
  // Initialize team topology and fan-in/out relationships:
  ThreadsExec::set_threads_relationships( team_topo , s_threads_exec );

  // Initial allocations:
  ThreadsExec::resize_reduce_scratch( 4096 );
  ThreadsExec::resize_shared_scratch( 4096 );
}

//----------------------------------------------------------------------------

void ThreadsExec::finalize()
{
  verify_is_process("ThreadsExec::finalize");

  fence();

  resize_reduce_scratch(0);
  resize_shared_scratch(0);

  const unsigned begin = s_threads_process.m_thread_size ? 1 : 0 ;

  for ( unsigned i = s_threads_count ; begin < i-- ; ) {

    if ( s_threads_exec[i] ) {

      s_threads_exec[i]->m_state = ThreadsExec::Terminating ;

      wait_yield( s_threads_process.m_state , ThreadsExec::Inactive );

      s_threads_process.m_state = ThreadsExec::Inactive ;
    }
  }

  if ( s_threads_process.m_thread_size ) {
    ( & s_threads_process )->~ThreadsExec();
    s_threads_exec[0] = 0 ;
    KokkosArray::hwloc::unbind_this_thread();
  }

  s_threads_count = 0 ;
}

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace Kokkos */



/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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

#ifndef KOKKOS_THREADSEXEC_HPP
#define KOKKOS_THREADSEXEC_HPP

#include <stdio.h>

#include <utility>
#include <impl/Kokkos_spinwait.hpp>

#include <Kokkos_Atomic.hpp>

#include <sys/types.h>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class > struct ThreadsExecAdapter ;

//----------------------------------------------------------------------------

class ThreadsExecTeamMember ;

class ThreadsExec {
public:

  // Fan array has log_2(NT) reduction threads plus 2 scan threads
  // Currently limited to 16k threads.
  enum { MAX_FAN_COUNT    = 16 };
  enum { MAX_THREAD_COUNT = 1 << ( MAX_FAN_COUNT - 2 ) };
  enum { VECTOR_LENGTH    = 8 };

  /** \brief States of a worker thread */
  enum { Terminating ///<  Termination in progress
       , Inactive    ///<  Exists, waiting for work
       , Active      ///<  Exists, performing work
       , Rendezvous  ///<  Exists, waiting in a barrier or reduce

       , ScanCompleted
       , ScanAvailable
       , ReductionAvailable
       };

private:

  friend class ThreadsExecTeamMember ;
  friend class Kokkos::Threads ;

  // Fan-in operations' root is the highest ranking thread
  // to place the 'scan' reduction intermediate values on
  // the threads that need them.
  // For a simple reduction the thread location is arbitrary.

  /** \brief  Reduction memory reserved for team reductions */
  enum { REDUCE_TEAM_BASE = 512 };

  ThreadsExec * const * m_pool_base ; ///< Base for pool fan-in
  ThreadsExec * const * m_team_base ; ///< Base for team fan-in

  pthread_t     m_pthread_id ;       ///< Pthread ID

  void        * m_alloc_reduce ;     ///< Reduction allocated memory
  void        * m_alloc_shared ;     ///< Team-shared allocated memory
  void        * m_team_shared ;      ///< Team-shared memory

  int           m_team_shared_end ;  ///< End of team-shared memory
  int           m_team_shared_iter ; ///< Current offset for team-shared memory

  int           m_pool_rank ;
  int           m_pool_size ;
  int           m_pool_fan_size ;

  int           m_team_rank ;
  int           m_team_size ;
  int           m_team_fan_size ;

  int           m_league_rank ;
  int           m_league_end ;
  int           m_league_size ;

  int volatile  m_pool_state ;  ///< State for global synchronizations
  int volatile  m_team_state ;  ///< State for team synchronizations

  static void global_lock();
  static void global_unlock();
  static bool spawn();

  static void execute_sleep( ThreadsExec & , const void * );
  static void execute_reduce_resize( ThreadsExec & , const void * );
  static void execute_shared_resize( ThreadsExec & , const void * );
  static void execute_get_binding(   ThreadsExec & , const void * );

  ThreadsExec( const ThreadsExec & );
  ThreadsExec & operator = ( const ThreadsExec & );

  static void execute_serial( void (*)( ThreadsExec & , const void * ) );

  inline void * reduce_team() const { return m_alloc_reduce ; }

  template < typename T >
  inline volatile T * team_reduce_value() const
    { return (volatile T *) m_alloc_reduce ; }

public:

  KOKKOS_INLINE_FUNCTION int pool_size() const { return m_pool_size ; }
  KOKKOS_INLINE_FUNCTION int pool_rank() const { return m_pool_rank ; }

  static int team_alloc( int team_size );

  static int get_thread_count();
  static ThreadsExec * get_thread( const int init_thread_rank );

  inline void * reduce_base() const { return ((unsigned char *) m_alloc_reduce) + REDUCE_TEAM_BASE ; }

  static void driver(void);

  void set_team_relations();

  ~ThreadsExec();
  ThreadsExec();

  static void resize_reduce_scratch( size_t );
  static void resize_shared_scratch( size_t );

  static void * root_reduce_scratch();

  static bool is_process();

  static void verify_is_process( const std::string & , const bool initialized );

  static int is_initialized();

  static void initialize( unsigned thread_count ,
                          unsigned use_numa_count ,
                          unsigned use_cores_per_numa ,
                          bool allow_asynchronous_threadpool );

  static void finalize();

  /* Given a requested team size, return valid team size */
  static unsigned team_size_valid( unsigned );

  static void print_configuration( std::ostream & , const bool detail = false );

  //------------------------------------

  static void wait_yield( volatile int & , const int );

  //------------------------------------
  // All-thread functions:

  inline
  std::pair< size_t , size_t >
  work_range( const size_t work_count ) const
  {
    typedef integral_constant< size_t , VECTOR_LENGTH - 1 > work_mask ;

    // work per thread rounded up and aligned to vector length:

    const size_t work_per_thread =
      ( ( ( work_count + m_pool_size - 1 ) / m_pool_size ) + work_mask::value ) & ~(work_mask::value);

    const size_t work_begin = std::min( work_count , work_per_thread * m_pool_rank );
    const size_t work_end   = std::min( work_count , work_per_thread + work_begin );

    return std::pair< size_t , size_t >( work_begin , work_end );
  }

  template< class Functor >
  inline
  void fan_in_reduce( const Functor & f ) const
    {
      typedef ReduceAdapter< Functor > Reduce ;

      const int rev_rank  = m_pool_size - ( m_pool_rank + 1 );

      for ( int i = 0 ; i < m_pool_fan_size ; ++i ) {

        ThreadsExec & fan = *m_pool_base[ rev_rank + ( 1 << i ) ] ;

        Impl::spinwait( fan.m_pool_state , ThreadsExec::Active );

        Reduce::join( f , reduce_base() , fan.reduce_base() );
      }

      if ( ! rev_rank ) {
        Reduce::final( f , reduce_base() );
      }
    }

  inline
  void fan_in() const
    {
      const int rev_rank = m_pool_size - ( m_pool_rank + 1 );

      for ( int i = 0 ; i < m_pool_fan_size ; ++i ) {
        Impl::spinwait( m_pool_base[rev_rank+(1<<i)]->m_pool_state , ThreadsExec::Active );
      }
    }

  template< class FunctorType >
  inline
  void scan_large( const FunctorType & f )

    {
      // Sequence of states:
      //  0) Active             : entry and exit state
      //  1) ReductionAvailable : reduction value available
      //  2) ScanAvailable      : inclusive scan value available
      //  3) Rendezvous         : All threads inclusive scan value are available
      //  4) ScanCompleted      : exclusive scan value copied

      typedef ReduceAdapter< FunctorType > Reduce ;
      typedef typename Reduce::scalar_type scalar_type ;

      const int      rev_rank = m_pool_size - ( m_pool_rank + 1 );
      const unsigned count    = Reduce::value_count( f );

      scalar_type * const work_value = (scalar_type *) reduce_base();

      //--------------------------------
      // Fan-in reduction with highest ranking thread as the root
      for ( int i = 0 ; i < m_pool_fan_size ; ++i ) {
        ThreadsExec & fan = *m_pool_base[ rev_rank + (1<<i) ];

        // Wait: Active -> ReductionAvailable (or ScanAvailable)
        Impl::spinwait( fan.m_pool_state , ThreadsExec::Active );
        Reduce::join( f , work_value , fan.reduce_base() );
      }

      // Copy reduction value to scan value before releasing from this phase.
      for ( unsigned i = 0 ; i < count ; ++i ) { work_value[i+count] = work_value[i] ; }

      if ( rev_rank ) {

        // Set: Active -> ReductionAvailable
        m_pool_state = ThreadsExec::ReductionAvailable ;

        // Wait for contributing threads' scan value to be available.
        if ( ( 1 << m_pool_fan_size ) < ( m_pool_rank + 1 ) ) {
          ThreadsExec & th = *m_pool_base[ rev_rank + ( 1 << m_pool_fan_size ) ] ;

          // Wait: Active             -> ReductionAvailable
          // Wait: ReductionAvailable -> ScanAvailable
          Impl::spinwait( th.m_pool_state , ThreadsExec::Active );
          Impl::spinwait( th.m_pool_state , ThreadsExec::ReductionAvailable );

          Reduce::join( f , work_value + count , ((scalar_type *)th.reduce_base()) + count );
        }

        // This thread has completed inclusive scan
        // Set: ReductionAvailable -> ScanAvailable
        m_pool_state = ThreadsExec::ScanAvailable ;

        // Wait for all threads to complete inclusive scan
        // Wait: ScanAvailable -> Rendezvous
        Impl::spinwait( m_pool_state , ThreadsExec::ScanAvailable );
      }

      //--------------------------------

      for ( int i = 0 ; i < m_pool_fan_size ; ++i ) {
        ThreadsExec & fan = *m_pool_base[ rev_rank + (1<<i) ];
        // Wait: ReductionAvailable -> ScanAvailable
        Impl::spinwait( fan.m_pool_state , ThreadsExec::ReductionAvailable );
        // Set: ScanAvailable -> Rendezvous
        fan.m_pool_state = ThreadsExec::Rendezvous ;
      }

      // All threads have completed the inclusive scan.
      // All non-root threads are in the Rendezvous state.
      // Threads are free to overwrite their reduction value.
      //--------------------------------

      if ( ( rev_rank + 1 ) < m_pool_size ) {
        // Exclusive scan: copy the previous thread's inclusive scan value

        ThreadsExec & th = *m_pool_base[ rev_rank + 1 ] ; // Not the root thread

        const scalar_type * const src_value = ((scalar_type *)th.reduce_base()) + count ;

        for ( unsigned j = 0 ; j < count ; ++j ) { work_value[j] = src_value[j]; }
      }
      else {
        (void) Reduce::init( f , work_value );
      }

      //--------------------------------
      // Wait for all threads to copy previous thread's inclusive scan value
      // Wait for all threads: Rendezvous -> ScanCompleted
      for ( int i = 0 ; i < m_pool_fan_size ; ++i ) {
        Impl::spinwait( m_pool_base[ rev_rank + (1<<i) ]->m_pool_state , ThreadsExec::Rendezvous );
      }
      if ( rev_rank ) {
        // Set: ScanAvailable -> ScanCompleted
        m_pool_state = ThreadsExec::ScanCompleted ;
        // Wait: ScanCompleted -> Active
        Impl::spinwait( m_pool_state , ThreadsExec::ScanCompleted );
      }
      // Set: ScanCompleted -> Active
      for ( int i = 0 ; i < m_pool_fan_size ; ++i ) {
        m_pool_base[ rev_rank + (1<<i) ]->m_pool_state = ThreadsExec::Active ;
      }
    }

  template< class FunctorType >
  inline
  void scan_small( const FunctorType & f )
    {
      typedef ReduceAdapter< FunctorType > Reduce ;
      typedef typename Reduce::scalar_type scalar_type ;

      const int      rev_rank = m_pool_size - ( m_pool_rank + 1 );
      const unsigned count    = Reduce::value_count( f );

      scalar_type * const work_value = (scalar_type *) reduce_base();

      //--------------------------------
      // Fan-in reduction with highest ranking thread as the root
      for ( int i = 0 ; i < m_pool_fan_size ; ++i ) {
        // Wait: Active -> Rendezvous
        Impl::spinwait( m_pool_base[ rev_rank + (1<<i) ]->m_pool_state , ThreadsExec::Active );
      }

      for ( unsigned i = 0 ; i < count ; ++i ) { work_value[i+count] = work_value[i]; }

      if ( rev_rank ) {
        m_pool_state = ThreadsExec::Rendezvous ;
        // Wait: Rendezvous -> Active
        Impl::spinwait( m_pool_state , ThreadsExec::Rendezvous );
      }
      else {
        // Root thread does the thread-scan before releasing threads

        scalar_type * ptr_prev = 0 ;

        for ( int rank = 0 ; rank < m_pool_size ; ++rank ) {
          scalar_type * const ptr = (scalar_type *) get_thread( rank )->reduce_base();
          if ( rank ) {
            for ( unsigned i = 0 ; i < count ; ++i ) { ptr[i] = ptr_prev[ i + count ]; }
            Reduce::join( f , ptr + count , ptr );
          }
          else {
            (void) Reduce::init( f , ptr );
          }
          ptr_prev = ptr ;
        }
      }

      for ( int i = 0 ; i < m_pool_fan_size ; ++i ) {
        m_pool_base[ rev_rank + (1<<i) ]->m_pool_state = ThreadsExec::Active ;
      }
    }

  //------------------------------------
  // Team-only functions:

  void * get_shmem( const int size );

  KOKKOS_INLINE_FUNCTION void team_barrier()
    {
      const int rev_rank = m_team_size - ( m_team_rank + 1 );

      for ( int i = 0 ; i < m_team_fan_size ; ++i ) {
        Impl::spinwait( m_team_base[ rev_rank + (1<<i) ]->m_pool_state , ThreadsExec::Active );
      }
      if ( rev_rank ) {
        m_pool_state = Rendezvous ;
        Impl::spinwait( m_pool_state , ThreadsExec::Rendezvous );
      }
      for ( int i = 0 ; i < m_team_fan_size ; ++i ) {
        m_team_base[ rev_rank + (1<<i) ]->m_pool_state = ThreadsExec::Active ;
      }
    }

  template< class ArgType >
  KOKKOS_INLINE_FUNCTION
  ArgType team_scan( const ArgType & value , ArgType * const global_accum = 0 )
    {
      // Sequence of m_team_state states:
      //  0) Inactive            : entry and exit state
      //  1) ReductionAvailable  : reduction value available, waiting for scan value
      //  2) ScanAvailable       : reduction value available, scan value available
      //  3) Rendezvous          : broadcasting global inter-team accumulation value

      // Make sure there is enough scratch space:
      typedef typename if_c< 2 * sizeof(ArgType) < REDUCE_TEAM_BASE , ArgType , void >::type type ;

      const int rev_rank = m_team_size - ( m_team_rank + 1 );

      type * const work_value = (type*) reduce_team();

      // ThreadsExec::Inactive == m_team_state

      work_value[0] = value ;
      memory_fence();

      // Fan-in reduction, wait for source thread to complete it's fan-in reduction.
      for ( int i = 0 ; i < m_team_fan_size ; ++i ) {
        ThreadsExec & th = *m_team_base[ rev_rank + (1<<i) ];

        // Wait for source thread to exit Inactive state.
        Impl::spinwait( th.m_team_state , ThreadsExec::Inactive );
        // Source thread is 'ReductionAvailable' or 'ScanAvailable'
        work_value[0] += ((volatile type*)th.reduce_team())[0];
        memory_fence();
      }

      work_value[1] = work_value[0] ;
      memory_fence();

      if ( rev_rank ) {

        m_team_state = ThreadsExec::ReductionAvailable ; // Reduction value is available.

        // Wait for contributing threads' scan value to be available.
        if ( ( 1 << m_team_fan_size ) < ( m_team_rank + 1 ) ) {
          ThreadsExec & th = *m_team_base[ rev_rank + ( 1 << m_team_fan_size ) ];

          // Wait: Inactive -> ReductionAvailable
          Impl::spinwait( th.m_team_state , ThreadsExec::Inactive );
          // Wait: ReductionAvailable -> ScanAvailable:
          Impl::spinwait( th.m_team_state , ThreadsExec::ReductionAvailable );

          work_value[1] += ((volatile type*)th.reduce_team())[1] ;
          memory_fence();
        }

        m_team_state = ThreadsExec::ScanAvailable ; // Scan value is available.
      }
      else {
         // Root thread add team's total to global inter-team accumulation
        work_value[0] = global_accum ? atomic_fetch_add( global_accum , work_value[0] ) : 0 ;
      }

      for ( int i = 0 ; i < m_team_fan_size ; ++i ) {
        ThreadsExec & th = *m_team_base[ rev_rank + (1<<i) ];
        // Wait: ReductionAvailable -> ScanAvailable
        Impl::spinwait( th.m_team_state , ThreadsExec::ReductionAvailable );
        // Wait: ScanAvailable -> Rendezvous
        Impl::spinwait( th.m_team_state , ThreadsExec::ScanAvailable );
      }

      // All fan-in threads are in the ScanAvailable state
      if ( rev_rank ) {
        m_team_state = ThreadsExec::Rendezvous ;
        Impl::spinwait( m_team_state , ThreadsExec::Rendezvous );
      }

      // Broadcast global inter-team accumulation value
      volatile type & global_val = work_value[0] ;
      for ( int i = 0 ; i < m_team_fan_size ; ++i ) {
        ThreadsExec & th = *m_team_base[ rev_rank + (1<<i) ];
        ((volatile type*)th.reduce_team())[0] = global_val ;
        memory_fence();
        th.m_team_state = ThreadsExec::Inactive ;
      }
      // Exclusive scan, subtract contributed value
      return global_val + work_value[1] - value ;
    }

  /*  When a functor using the 'device' interface requests
   *  more teams than are initialized the parallel operation
   *  must loop over a range of league ranks with a team_barrier
   *  between each iteration.
   */
  bool team_work_avail()
    { m_team_shared_iter = 0 ; return m_league_rank < m_league_end ; }

  void team_work_next()
    { if ( ++m_league_rank < m_league_end ) team_barrier(); }

  //------------------------------------
  /** \brief  Wait for previous asynchronous functor to
   *          complete and release the Threads device.
   *          Acquire the Threads device and start this functor.
   */
  static void start( void (*)( ThreadsExec & , const void * ) , const void * ,
                     int work_league_size = 0 ,
                     int work_team_size = 0 );

  static unsigned team_max();
  static unsigned team_recommended();
  static unsigned hardware_thread_id();
  static unsigned max_hardware_threads();

  static int  in_parallel();
  static void fence();
  static bool sleep();
  static bool wake();
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class ThreadsExecTeamMember {
private:

  Impl::ThreadsExec & m_exec ;

  
  
  typedef Kokkos::Threads execution_space ;

  // Fan-in and wait until the matching fan-out is called.
  // The root thread which does not wait will return true.
  // All other threads will return false during the fan-out.
  KOKKOS_INLINE_FUNCTION bool team_fanin() const
    {
      const int rev_rank = m_exec.m_team_size - ( m_exec.m_team_rank + 1 );
      const bool is_root = ! rev_rank ;

      int n , j ;

      // Wait for fan-in threads
      for ( n = 1 ; ( ! ( rev_rank & n ) ) && ( ( j = rev_rank + n ) < m_exec.m_team_size ) ; n <<= 1 ) {
        Impl::spinwait( m_exec.m_team_base[j]->m_pool_state , ThreadsExec::Active );
      }

      // If not root then wait for release
      if ( ! is_root ) {
        m_exec.m_pool_state = ThreadsExec::Rendezvous ;
        Impl::spinwait( m_exec.m_pool_state , ThreadsExec::Rendezvous );
      }

      return is_root ;
    }

  KOKKOS_INLINE_FUNCTION void team_fanout() const
    {
      const int rev_rank = m_exec.m_team_size - ( m_exec.m_team_rank + 1 );
      int n , j ;
      for ( n = 1 ; ( ! ( rev_rank & n ) ) && ( ( j = rev_rank + n ) < m_exec.m_team_size ) ; n <<= 1 ) {
        m_exec.m_team_base[j]->m_pool_state = ThreadsExec::Active ;
      }
    }

public:

  KOKKOS_INLINE_FUNCTION
  execution_space::scratch_memory_space  team_shmem() const
    { return execution_space::scratch_memory_space( m_exec ); }

  KOKKOS_INLINE_FUNCTION int league_rank() const { return m_exec.m_league_rank ; }
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_exec.m_league_size ; }
  KOKKOS_INLINE_FUNCTION int team_rank() const { return m_exec.m_team_rank ; }
  KOKKOS_INLINE_FUNCTION int team_size() const { return m_exec.m_team_size ; }

  KOKKOS_INLINE_FUNCTION void team_barrier() const
    {
      team_fanin();
      team_fanout();
    }

  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_reduce( const Type & value ) const
    {
      // Make sure there is enough scratch space:
      typedef typename if_c< sizeof(Type) < ThreadsExec::REDUCE_TEAM_BASE , Type , void >::type type ;

      *((volatile type*) m_exec.reduce_team() ) = value ;

      memory_fence();

      type & accum = *((type *) m_exec.m_team_base[0]->reduce_team() );

      if ( team_fanin() ) {
        for ( int i = 1 ; i < m_exec.m_team_size ; ++i ) {
          accum += *((type *) m_exec.m_team_base[i]->reduce_team() );
        }
        memory_fence();
      }

      team_fanout();

      return accum ;
    }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering
   *          with intra-team non-deterministic ordering accumulation.
   *
   *  The global inter-team accumulation value will, at the end of the
   *  league's parallel execution, be the scan's total.
   *  Parallel execution ordering of the league's teams is non-deterministic.
   *  As such the base value for each team's scan operation is similarly
   *  non-deterministic.
   */
  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value , Type * const global_accum ) const
    {
      // Make sure there is enough scratch space:
      typedef typename if_c< sizeof(Type) < ThreadsExec::REDUCE_TEAM_BASE , Type , void >::type type ;

      *((volatile type*) m_exec.reduce_team() ) = value ;

      memory_fence();

      if ( team_fanin() ) {
        type accum = *((type *) m_exec.m_team_base[0]->reduce_team() );

        // Copy from lower to higher rank team member: { N-1 , N-2 , N-3 , ... , 0 }
        for ( int i = 1 ; i < m_exec.m_team_size ; ++i ) {
          accum += ( *((type *) m_exec.m_team_base[i-1]->reduce_team() ) =
                     *((type *) m_exec.m_team_base[i  ]->reduce_team() ) );
        }

        *((type *) m_exec.m_team_base[ m_exec.m_team_size - 1 ]->reduce_team() ) =
          global_accum ? atomic_fetch_add( global_accum , accum ) : 0 ;

        // Join from lower rank to higher rank

        for ( int i = m_exec.m_team_size ; --i ; ) {
          *((type *) m_exec.m_team_base[i-1]->reduce_team() ) +=
          *((type *) m_exec.m_team_base[i  ]->reduce_team() );
        }

        memory_fence();
      }

      team_fanout();

      return *((volatile type*) m_exec.reduce_team() );
    }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value ) const
    { return m_exec.template team_scan<Type>( value , 0 ); }

  //----------------------------------------
  // Private for the driver

  template< class WorkArgTag >
  ThreadsExecTeamMember( Impl::ThreadsExec & exec , const TeamPolicy< execution_space , WorkArgTag > & team )
    : m_exec( exec )
    {}

  void reset_scratch_space()
    { m_exec.m_team_shared_iter = 0 ; }

  bool valid_team() const
    { return m_exec.m_league_rank < m_exec.m_league_end ; }

  void next_team()
    { if ( ++m_exec.m_league_rank < m_exec.m_league_end ) m_exec.team_barrier(); }
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

inline int Threads::in_parallel()
{ return Impl::ThreadsExec::in_parallel(); }

inline int Threads::is_initialized()
{ return Impl::ThreadsExec::is_initialized(); }

inline void Threads::initialize(
  unsigned threads_count ,
  unsigned use_numa_count ,
  unsigned use_cores_per_numa ,
  bool allow_asynchronous_threadpool )
{
  Impl::ThreadsExec::initialize( threads_count , use_numa_count , use_cores_per_numa , allow_asynchronous_threadpool );
}

inline void Threads::finalize()
{
  Impl::ThreadsExec::finalize();
}

inline void Threads::print_configuration( std::ostream & s , const bool detail )
{
  Impl::ThreadsExec::print_configuration( s , detail );
}

KOKKOS_INLINE_FUNCTION unsigned Threads::team_max()
{ return Impl::ThreadsExec::team_max() ; }

KOKKOS_INLINE_FUNCTION unsigned Threads::team_recommended()
{ return Impl::ThreadsExec::team_recommended() ; }

KOKKOS_INLINE_FUNCTION unsigned Threads::hardware_thread_id()
{ return Impl::ThreadsExec::hardware_thread_id() ; }

KOKKOS_INLINE_FUNCTION unsigned Threads::max_hardware_threads()
{ return Impl::ThreadsExec::max_hardware_threads() ; }

inline bool Threads::sleep()
{ return Impl::ThreadsExec::sleep() ; }

inline bool Threads::wake()
{ return Impl::ThreadsExec::wake() ; }

inline void Threads::fence()
{ Impl::ThreadsExec::fence() ; }

KOKKOS_INLINE_FUNCTION int Threads::league_rank() const
{ return m_exec.m_league_rank ; }

KOKKOS_INLINE_FUNCTION int Threads::league_size() const
{ return m_exec.m_league_size ; }

KOKKOS_INLINE_FUNCTION int Threads::team_rank() const
{ return m_exec.m_team_rank ; }

KOKKOS_INLINE_FUNCTION int Threads::team_size() const
{ return m_exec.m_team_size ; }

KOKKOS_INLINE_FUNCTION void Threads::team_barrier()
{ return m_exec.team_barrier(); }

inline Threads::Threads( Impl::ThreadsExec & t ) : m_exec( t ) {}

template< typename Type >
KOKKOS_INLINE_FUNCTION Type Threads::team_scan( const Type & value )
{ return m_exec.team_scan( value ); }

template< typename TypeLocal , typename TypeGlobal >
KOKKOS_INLINE_FUNCTION TypeGlobal Threads::team_scan( const TypeLocal & value , TypeGlobal * const global_accum )
{ return m_exec.template team_scan< TypeGlobal >( value , global_accum ); }

KOKKOS_INLINE_FUNCTION
void * Threads::get_shmem( const int size ) const { return m_exec.get_shmem( size ); }

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template < class WorkArgTag >
class TeamPolicy< Kokkos::Threads , WorkArgTag > {
private:

  const int m_league_size ;
  const int m_team_size ;
  const int m_team_alloc ;

public:

  typedef Impl::ExecutionPolicyTag   kokkos_tag ;      ///< Concept tag
  typedef Kokkos::Threads            execution_space ; ///< Execution space

  inline int team_size() const { return m_team_size ; }
  inline int league_size() const { return m_league_size ; }

  /** \brief  Specify league size, request team size */
  TeamPolicy( execution_space & , int league_size_request , int team_size_request )
    : m_league_size( league_size_request )
    , m_team_size( team_size_request < int(execution_space::team_max())
                 ? team_size_request : int(execution_space::team_max()) )
    , m_team_alloc( Impl::ThreadsExec::team_alloc( m_team_size ) )
    { }

  TeamPolicy( int league_size_request , int team_size_request )
    : m_league_size( league_size_request )
    , m_team_size( team_size_request < int(execution_space::team_max())
                 ? team_size_request : int(execution_space::team_max()) )
    , m_team_alloc( Impl::ThreadsExec::team_alloc( m_team_size ) )
    { }

  typedef Impl::ThreadsExecTeamMember member_type ;

  friend class Impl::ThreadsExecTeamMember ;
};

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #define KOKKOS_THREADSEXEC_HPP */


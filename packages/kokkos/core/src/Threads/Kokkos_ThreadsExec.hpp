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

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class > struct ThreadsExecAdapter ;

//----------------------------------------------------------------------------

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

       , ScanAvailable
       , ReductionAvailable
       };

private:

  friend class Kokkos::Threads ;

  // Fan-in operations' root is the highest ranking thread
  // to place the 'scan' reduction intermediate values on
  // the threads that need them.
  // For a simple reduction the thread location is arbitrary.

  /** \brief  Reduction memory reserved for team reductions */
  enum { REDUCE_TEAM_BASE = 512 };

  void        * m_reduce ;      ///< Reduction memory
  void        * m_shared ;      ///< Team-shared memory
  int           m_shared_end ;  ///< End of team-shared memory
  int           m_shared_iter ; ///< Current offset for team-shared memory
  int volatile  m_state ;       ///< State for global synchronizations
  int volatile  m_state_team ;  ///< State for team synchronizations
  int           m_fan_size ;
  int           m_fan_team_size ;

  int           m_team_rank ;
  int           m_team_size ;
  int           m_init_league_rank ;
  int           m_init_league_size ;
  int           m_init_thread_rank ;
  int           m_init_thread_size ;

  int           m_work_league_rank ;
  int           m_work_league_end ;
  int           m_work_league_size ;

  ThreadsExec * m_fan[ MAX_FAN_COUNT ] ;
  ThreadsExec * m_fan_team[ MAX_FAN_COUNT ] ;

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

  inline void * reduce_team() const { return m_reduce ; }

public:

  static int get_thread_count();
  static ThreadsExec * get_thread( const int init_thread_rank );

  inline void * reduce_base() const { return ((unsigned char *) m_reduce) + REDUCE_TEAM_BASE ; }

  static void driver(void);

  ~ThreadsExec();
  ThreadsExec();

  static void set_threads_relationships( const std::pair<unsigned,unsigned> team_topo );

  static void resize_reduce_scratch( size_t );
  static void resize_shared_scratch( size_t );

  static void * root_reduce_scratch();

  static bool is_process();

  static void verify_is_process( const std::string & , const bool initialized );

  static int is_initialized();

  static void initialize( const std::pair<unsigned,unsigned> team_topo ,
                                std::pair<unsigned,unsigned> core_topo );

  static void finalize();

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
      ( ( ( work_count + m_init_thread_size - 1 ) / m_init_thread_size ) + work_mask::value ) & ~(work_mask::value);

    const size_t work_begin = std::min( work_count , work_per_thread * m_init_thread_rank );
    const size_t work_end   = std::min( work_count , work_per_thread + work_begin );

    return std::pair< size_t , size_t >( work_begin , work_end );
  }

  template< class Functor >
  inline
  void fan_in_reduce( const Functor & f ) const
    {
      typedef ReduceAdapter< Functor > Reduce ;

      for ( int i = 0 ; i < m_fan_size ; ++i ) {

        ThreadsExec & fan = *m_fan[i] ;

        Impl::spinwait( fan.m_state , ThreadsExec::Active );

        f.join( Reduce::reference( reduce_base() ) ,
                Reduce::reference( fan.reduce_base() ) );
      }
    }

  inline
  void fan_in() const
    {
      for ( int i = 0 ; i < m_fan_size ; ++i ) {
        Impl::spinwait( m_fan[i]->m_state , ThreadsExec::Active );
      }
    }

  template< class FunctorType >
  inline
  void scan_large( const FunctorType & f )
    {
      // Sequence of states:
      //  0) Active             : entry and exit state
      //  1) ReductionAvailable : reduction value available
      //  1) Rendezvous         : all reduction values available and copied
      //  2) ScanAvailable      : scan value available

      typedef ReduceAdapter< FunctorType > Reduce ;
      typedef typename Reduce::scalar_type scalar_type ;

      const bool     not_root = m_init_thread_rank + 1 < m_init_thread_size ;
      const unsigned count    = Reduce::value_count( f );

      scalar_type * const work_value = (scalar_type *) reduce_base();

      //--------------------------------
      // Fan-in reduction with highest ranking thread as the root
      for ( int i = 0 ; i < m_fan_size ; ++i ) {
        ThreadsExec & fan = *m_fan[i];

        // Wait: Active -> ReductionAvailable
        Impl::spinwait( fan.m_state , ThreadsExec::Active );
        f.join( Reduce::reference( work_value ) , Reduce::reference( fan.reduce_base() ) );
      }

      // Copy reduction value to scan value before releasing from this phase.
      for ( unsigned i = 0 ; i < count ; ++i ) { work_value[i+count] = work_value[i] ; }

      if ( not_root ) {
        m_state = ThreadsExec::ReductionAvailable ;
        // Wait: ReductionAvailable -> Rendezvous
        Impl::spinwait( m_state , ThreadsExec::ReductionAvailable );
      }

      for ( int i = 0 ; i < m_fan_size ; ++i ) {
        m_fan[i]->m_state = ThreadsExec::Rendezvous ;
      }

      // All non-root threads are now in the Rendezvous state
      //--------------------------------

      if ( not_root ) {

        // Wait for contributing threads' scan value to be available.
        if ( m_fan[ m_fan_size ] ) {
          ThreadsExec & th = *m_fan[ m_fan_size ] ;

          // Wait: Rendezvous -> ScanAvailable
          Impl::spinwait( th.m_state , ThreadsExec::Rendezvous );

          f.join( Reduce::reference( work_value + count ) ,
                  Reduce::reference( ((scalar_type *)th.reduce_base()) + count ) );
        }

        m_state = ThreadsExec::ScanAvailable ;
      }

      //--------------------------------

      if ( m_fan[ m_fan_size + 1 ] ) {
        ThreadsExec & th = *m_fan[ m_fan_size + 1 ] ; // Not the root thread

        // Wait: Rendezvous -> ScanAvailable
        Impl::spinwait( th.m_state , ThreadsExec::Rendezvous );

        const scalar_type * const src_value = ((scalar_type *)th.reduce_base()) + count ;

        for ( unsigned j = 0 ; j < count ; ++j ) { work_value[j] = src_value[j]; }

        th.m_state = ThreadsExec::Active ; // Release the source thread
      }
      else {
        f.init( Reduce::reference( work_value ) );
      }

      // Wait for scan value to be claimed before exiting.
      Impl::spinwait( m_state , ThreadsExec::ScanAvailable );
    }

  template< class FunctorType >
  inline
  void scan_small( const FunctorType & f )
    {
      typedef ReduceAdapter< FunctorType > Reduce ;
      typedef typename Reduce::scalar_type scalar_type ;

      const bool     not_root = m_init_thread_rank + 1 < m_init_thread_size ;
      const unsigned count    = Reduce::value_count( f );

      scalar_type * const work_value = (scalar_type *) reduce_base();

      //--------------------------------
      // Fan-in reduction with highest ranking thread as the root
      for ( int i = 0 ; i < m_fan_size ; ++i ) {
        // Wait: Active -> Rendezvous
        Impl::spinwait( m_fan[i]->m_state , ThreadsExec::Active );
      }

      for ( unsigned i = 0 ; i < count ; ++i ) { work_value[i+count] = work_value[i]; }

      if ( not_root ) {
        m_state = ThreadsExec::Rendezvous ;
        // Wait: Rendezvous -> Active
        Impl::spinwait( m_state , ThreadsExec::Rendezvous );
      }
      else {
        // Root thread does the thread-scan before releasing threads

        scalar_type * ptr_prev = 0 ;

        for ( int rank = 0 ; rank < m_init_thread_size ; ++rank ) {
          scalar_type * const ptr = (scalar_type *) get_thread( rank )->reduce_base();
          if ( rank ) {
            for ( unsigned i = 0 ; i < count ; ++i ) { ptr[i] = ptr_prev[ i + count ]; }
            f.join( Reduce::reference( ptr + count ), Reduce::reference( ptr ) );
          }
          else {
            f.init( Reduce::reference( ptr ) );
          }
          ptr_prev = ptr ;
        }
      }

      for ( int i = 0 ; i < m_fan_size ; ++i ) {
        m_fan[i]->m_state = ThreadsExec::Active ;
      }
    }

  //------------------------------------
  // Team-only functions:

  void * get_shmem( const int size );

  void team_barrier()
    {
      const bool not_root = m_team_rank + 1 < m_team_size ;

      for ( int i = 0 ; i < m_fan_team_size ; ++i ) {
        Impl::spinwait( m_fan_team[i]->m_state , ThreadsExec::Active );
      }
      if ( not_root ) {
        m_state = Rendezvous ;
        Impl::spinwait( m_state , ThreadsExec::Rendezvous );
      }
      for ( int i = 0 ; i < m_fan_team_size ; ++i ) {
        m_fan_team[i]->m_state = ThreadsExec::Active ;
      }
    }

  template< class ArgType >
  inline
  ArgType team_scan( const ArgType & value , ArgType * const global_accum = 0 )
    {
      // Sequence of m_state_team states:
      //  0) Inactive            : entry and exit state
      //  1) ReductionAvailable  : reduction value available, waiting for scan value
      //  2) ScanAvailable       : reduction value available, scan value available
      //  3) Rendezvous          : broadcasting global iinter-team accumulation value

      // Make sure there is enough scratch space:
      typedef typename if_c< 2 * sizeof(ArgType) < REDUCE_TEAM_BASE , ArgType , void >::type type ;

      const bool not_root = m_team_rank + 1 < m_team_size ;
      type * const work_value = (type*) reduce_team();

      // ThreadsExec::Inactive == m_state_team

      work_value[0] = value ;

      // Fan-in reduction, wait for source thread to complete it's fan-in reduction.
      for ( int i = 0 ; i < m_fan_team_size ; ++i ) {
        ThreadsExec & th = *m_fan_team[i];

        // Wait for source thread to exit Inactive state.
        Impl::spinwait( th.m_state_team , ThreadsExec::Inactive );
        // Source thread is 'ReductionAvailable' or 'ScanAvailable'
        work_value[0] += ((volatile type*)th.reduce_team())[0];
      }

      work_value[1] = work_value[0] ;

      if ( not_root ) {

        m_state_team = ThreadsExec::ReductionAvailable ; // Reduction value is available.

        // Wait for contributing threads' scan value to be available.
        if ( m_fan_team[ m_fan_team_size ] ) {
          ThreadsExec & th = *m_fan_team[ m_fan_team_size ] ;

          // Wait: Inactive -> ReductionAvailable
          Impl::spinwait( th.m_state_team , ThreadsExec::Inactive );
          // Wait: ReductionAvailable -> ScanAvailable:
          Impl::spinwait( th.m_state_team , ThreadsExec::ReductionAvailable );

          work_value[1] += ((volatile type*)th.reduce_team())[1] ;
        }

        m_state_team = ThreadsExec::ScanAvailable ; // Scan value is available.
      }
      else {
         // Root thread add team's total to global inter-team accumulation
        work_value[0] = global_accum ? atomic_fetch_add( global_accum , work_value[0] ) : 0 ;
      }

      for ( int i = 0 ; i < m_fan_team_size ; ++i ) {
        ThreadsExec & th = *m_fan_team[i];
        // Wait: ReductionAvailable -> ScanAvailable
        Impl::spinwait( th.m_state_team , ThreadsExec::ReductionAvailable );
        // Wait: ScanAvailable -> Rendezvous
        Impl::spinwait( th.m_state_team , ThreadsExec::ScanAvailable );
      }

      // All fan-in threads are in the ScanAvailable state
      if ( not_root ) {
        m_state_team = ThreadsExec::Rendezvous ;
        Impl::spinwait( m_state_team , ThreadsExec::Rendezvous );
      }

      // Broadcast global inter-team accumulation value
      volatile type & global_val = work_value[0] ;
      for ( int i = 0 ; i < m_fan_team_size ; ++i ) {
        ThreadsExec & th = *m_fan_team[i];
        ((volatile type*)th.reduce_team())[0] = global_val ;
        th.m_state_team = ThreadsExec::Inactive ;
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
    {
      m_shared_iter = 0 ;
      return m_work_league_rank < m_work_league_end ;
    }

  void team_work_next()
    { if ( ++m_work_league_rank < m_work_league_end ) team_barrier(); }

  //------------------------------------
  /** \brief  Wait for previous asynchronous functor to
   *          complete and release the Threads device.
   *          Acquire the Threads device and start this functor.
   */
  static void start( void (*)( ThreadsExec & , const void * ) , const void * , int = 0 );

  static int league_max();
  static int team_max();

  static int  in_parallel();
  static void fence();
  static bool sleep();
  static bool wake();
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
  unsigned team_count ,
  unsigned threads_per_team ,
  unsigned use_numa_count ,
  unsigned use_cores_per_numa )
{
  Impl::ThreadsExec::initialize(
    std::pair<unsigned,unsigned>( team_count , threads_per_team ),
    std::pair<unsigned,unsigned>( use_numa_count , use_cores_per_numa ) );
}

inline void Threads::finalize()
{
  Impl::ThreadsExec::finalize();
}

inline void Threads::print_configuration( std::ostream & s , const bool detail )
{
  Impl::ThreadsExec::print_configuration( s , detail );
}

inline unsigned Threads::league_max()
{ return Impl::ThreadsExec::league_max() ; }

inline unsigned Threads::team_max()
{ return Impl::ThreadsExec::team_max() ; }

inline bool Threads::sleep()
{ return Impl::ThreadsExec::sleep() ; }

inline bool Threads::wake()
{ return Impl::ThreadsExec::wake() ; }

inline void Threads::fence()
{ Impl::ThreadsExec::fence() ; }

inline int Threads::league_rank() const
{ return m_exec.m_work_league_rank ; }

inline int Threads::league_size() const
{ return m_exec.m_work_league_size ; }

inline int Threads::team_rank() const
{ return m_exec.m_team_rank ; }

inline int Threads::team_size() const
{ return m_exec.m_team_size ; }

inline void Threads::team_barrier()
{ return m_exec.team_barrier(); }

inline Threads::Threads( Impl::ThreadsExec & t ) : m_exec( t ) {}

template< typename Type >
inline Type Threads::team_scan( const Type & value )
{ return m_exec.team_scan( value ); }

template< typename TypeLocal , typename TypeGlobal >
inline TypeGlobal Threads::team_scan( const TypeLocal & value , TypeGlobal * const global_accum )
{ return m_exec.template team_scan< TypeGlobal >( value , global_accum ); }

inline
void * Threads::get_shmem( const int size ) { return m_exec.get_shmem( size ); }

} /* namespace Kokkos */

#endif /* #define KOKKOS_THREADSEXEC_HPP */


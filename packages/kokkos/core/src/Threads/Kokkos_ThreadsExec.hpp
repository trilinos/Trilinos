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

#include <utility>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class > struct ThreadsExecAdapter ;

//----------------------------------------------------------------------------

class ThreadsExec {
public:

  enum { MAX_FAN_COUNT    = 16 };
  enum { MAX_THREAD_COUNT = 1 << MAX_FAN_COUNT };

  /** \brief States of a worker thread */
  enum { Terminating ///<  Termination in progress
       , Inactive    ///<  Exists, waiting for work
       , Active      ///<  Exists, performing work
       , Rendezvous  ///<  Exists, waiting in a barrier or reduce
       };

private:

  friend class Kokkos::Threads ;

  void        * m_reduce ;    ///< Reduction memory
  void        * m_shared ;    ///< Shared memory
  int           m_shared_end ;
  int           m_shared_iter ;
  int volatile  m_state ;
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

public:

  static void driver(void);

  ~ThreadsExec();
  ThreadsExec();

  static void set_threads_relationships( const std::pair<unsigned,unsigned> team_topo ,
                                         ThreadsExec * threads[] );

  static void resize_reduce_scratch( size_t );
  static void resize_shared_scratch( size_t );

  static void * root_reduce_scratch();

  static bool is_process();

  static void verify_is_process( const std::string & , const bool initialized );

  static void initialize( const std::pair<unsigned,unsigned> team_topo ,
                                std::pair<unsigned,unsigned> core_topo );

  static void finalize();

  static void print_configuration( std::ostream & , const bool detail = false );

  //------------------------------------

  static void wait( volatile int & , const int );
  static void wait_yield( volatile int & , const int );

  void * get_shmem( const int size );

  template< class FunctorType >
  typename ReduceAdapter< FunctorType >::reference_type
  reduce_value( const FunctorType & ) const
    { return ReduceAdapter< FunctorType >::reference( m_reduce ); }

  template< class Functor >
  inline
  typename enable_if< FunctorHasJoin< Functor >::value >::type
  fan_in( const Functor & f ) const
    {
      typedef ReduceAdapter< Functor > Reduce ;

      for ( int i = 0 ; i < m_fan_size ; ++i ) {

        ThreadsExec & fan = *m_fan[i] ;

        ThreadsExec::wait( fan.m_state , ThreadsExec::Active );

        f.join( Reduce::reference( m_reduce ) ,
                Reduce::reference( fan.m_reduce ) );
      }
    }

  template< class Functor >
  inline
  typename enable_if< ! FunctorHasJoin< Functor >::value >::type
  fan_in( const Functor & ) const
    {
      for ( int i = 0 ; i < m_fan_size ; ++i ) {
        ThreadsExec::wait( m_fan[i]->m_state , ThreadsExec::Active );
      }
    }

  void team_barrier()
    {
      for ( int i = 0 ; i < m_fan_team_size ; ++i ) {
        ThreadsExec::wait( m_fan[i]->m_state , ThreadsExec::Active );
      }
      if ( m_team_rank ) {
        m_state = Rendezvous ;
        ThreadsExec::wait( m_state , ThreadsExec::Rendezvous );
      }
      for ( int i = 0 ; i < m_fan_team_size ; ++i ) {
        m_fan[i]->m_state = ThreadsExec::Active ;
      }
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

  inline
  std::pair< size_t , size_t >
  work_range( const size_t work_count ) const
  {
    enum { work_align = Kokkos::HostSpace::WORK_ALIGNMENT };
    enum { work_shift = Kokkos::Impl::power_of_two< work_align >::value };
    enum { work_mask  = work_align - 1 };

    // unit_of_work_count = ( work_count + work_mask ) >> work_shift
    // unit_of_work_per_thread = ( unit_of_work_count + thread_count - 1 ) / thread_count
    // work_per_thread = unit_of_work_per_thread * work_align

    const size_t work_per_thread =
      ( ( ( ( work_count + work_mask ) >> work_shift ) + m_init_thread_size - 1 ) / m_init_thread_size ) << work_shift ;

    const size_t work_begin = std::min( m_init_thread_rank * work_per_thread , work_count );
    const size_t work_end   = std::min( work_begin + work_per_thread , work_count );

    return std::pair< size_t , size_t >( work_begin , work_end );
  }


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

inline void Threads::initialize( 
  const std::pair<unsigned,unsigned> league_team ,
  const std::pair<unsigned,unsigned> hardware_topology )
{
  Impl::ThreadsExec::initialize( league_team , hardware_topology );
}

inline void Threads::finalize()
{
  Impl::ThreadsExec::finalize();
}

inline void Threads::print_configuration( std::ostream & s , bool detail )
{
  Impl::ThreadsExec::print_configuration( s , detail );
}

inline int Threads::league_max()
{ return Impl::ThreadsExec::league_max() ; }

inline int Threads::team_max()
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

template< class ViewType >
inline typename ViewType::value_type
Threads::unordered_scan
   (typename ViewType::value_type& value, ViewType& scratch_view)
{
  //gives back an ordered sum within a team, last thread has the highest value)
  typename ViewType::value_type sum = team_scan(value);

  //
  typename ViewType::value_type global_sum;

  if(team_rank()==team_size()-1)
    global_sum = atomic_fetch_add(&scratch_view(0),sum+value);

  typename ViewType::value_type* temp =
    (typename ViewType::value_type *) get_shmem( sizeof(typename ViewType::value_type) );

  if(team_rank()==team_size()-1)
    temp[0]=global_sum;

  //TODO: reduce shared memory counter again?
  team_barrier();
  return temp[0] + sum;
}

//TODO: find optimal implementation CRT i.e.
//team_scan calculates the exclusive prefix sum of the values of the threads
//return values are monotonous increasing with thread index
template< typename T >
inline T Threads::team_scan( T& value ) {

  if( team_size() == 1 ) return 0;

  T* temp = (T*) get_shmem( sizeof(T) * team_size() );
  temp[team_rank()] = value;


  team_barrier();

  // If team sizes are small it is faster to use a single thread to do the scan
  if( team_rank() == 0 ) {
    T old = 0;
    for( int i = 0; i < team_size(); i++ ) {
    const T next = temp[i];
    temp[i] = old;
    old += next;
    }
  }
  team_barrier();

  //TODO: reduce sharedmem counter again? CRT
  return temp[team_rank()];

  //TODO: implement team_scan with log(n) algorithm? Probably not usefull if team_size is small CRT
}

inline
void * Threads::get_shmem( const int size ) { return m_exec.get_shmem( size ); }

} /* namespace Kokkos */

#endif /* #define KOKKOS_THREADSEXEC_HPP */

